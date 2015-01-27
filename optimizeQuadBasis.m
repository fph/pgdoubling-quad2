function [quad,invcond,swaps]=optimizeQuadBasis(quad,diagonalThreshold,offDiagonalThreshold,maxSwaps)
% given a symBasis, reduces it so that all elements are below a threshold
%
% [quad,swaps,invcond]=optimizeQuadBasis(quad,diagonalThreshold,offdiagonalThreshold,maxSwaps)
%
% invcond is the inverse condition number of all the performed
% transformations composed. If it is larger than some moderate value,
% consider recomputing the basis
%
% diagonalThreshold and offDiagonalThreshold are hard limits on the
% possible magnitudes of the elements of symBasisFromQuadBasis(quad).
%
% TODO: the condition number probably isn't the right thing --- (see other
% bases)
% conditioned
%
% (c) 2015-2015 F. Poloni <poloni@math.tu-berlin.de> and others
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

if not(exist('diagonalThreshold','var')) || isempty(diagonalThreshold)
    diagonalThreshold=2;
end

if not(exist('offDiagonalThreshold','var')) || isempty(offDiagonalThreshold)
    offDiagonalThreshold=diagonalThreshold + 1;
end

if diagonalThreshold<1+sqrt(eps(class(quad.X)))
    error('PGDoubling:thresholdTooSmall','you can only hope to enforce thresholds S=1+sqrt(eps) or larger');
end

if offDiagonalThreshold<diagonalThreshold
    error('PGDoubling:thresholdTooSmall','you can only hope to enforce thresholds T=S or larger');
end

if not(exist('maxSwaps','var')) || isempty(maxSwaps)
    maxSwaps=20*size(quad.X,2);
end

%optimization
swaps=0;
invcond=1;
indices = 1:length(quad.X); %helper for Matlab indexing
while(true)
    % TODO: use Natasa's heuristic instead of recomputing norms each time
    % row norms of A
    squaredRowNorms = sum(abs(quad.X(quad.v,quad.v)).^2,1);
    [maxSquaredRowNorm,maxi] = max(squaredRowNorms);
    if maxSquaredRowNorm > diagonalThreshold
        % I can't find how to avoid this index juggling --federico
        indicesInV = indices(quad.v);
        maxi = indicesInV(maxi);
        [quad.X, quad.v, stepcond] = updateQuadBasisOut(quad.X, quad.v, maxi);
        swaps = swaps + 1;
    else
        squaredColumnNorms = sum(abs(quad.X(~quad.v,~quad.v)).^2,2);
        [maxColumnNorm,maxj] = max(squaredColumnNorms);
        if maxColumnNorm > diagonalThreshold
            indicesNotInV = indices(~quad.v);
            maxj = indicesNotInV(maxj);
            [quad.X, quad.v, stepcond] = updateQuadBasisIn(quad.X, quad.v, maxj);
            swaps = swaps + 1;
        else
            dets = abs(quad.X(~quad.v,quad.v)).^2 + squaredColumnNorms*squaredRowNorms;
            [maxvec, maxis] = max(dets);
            [maxval, maxj] = max(maxvec);
            if maxval > offDiagonalThreshold^2
                indicesNotInV = indices(~quad.v);
                maxi = indicesNotInV(maxis(maxj));
                indicesInV = indices(quad.v);
                maxj = indicesInV(maxj);
                [quad.X, quad.v, stepcond] = updateQuadBasisInOut(quad.X, quad.v, maxi, maxj);
                swaps = swaps + 2;
            else
                break %finished optimization
            end
        end
    end
    invcond = invcond*stepcond;
    if swaps >= maxSwaps
        warning('PGDoubling:stagnated','failed to produce a X with elements below the required threshold. Try running with a larger threshold.');
        break
    end
end
