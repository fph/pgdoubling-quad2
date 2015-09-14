function [absdets,maxnorms]=optimizeQuadBasis(quad,threshold,maxSwaps)
% given a symBasis, reduces it so that all elements are below a threshold
%
% A version of optimizeQuadBasis that returns the sequence of |det(X)| and
% max |X| along the steps, for plotting.
%
% (c) 2015-2015 F. Poloni <poloni@math.tu-berlin.de> and others
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

if not(exist('threshold','var')) || isempty(threshold)
    threshold=2;
end

if threshold<1+sqrt(eps(class(quad.X)))
    error('PGDoubling:thresholdTooSmall','you can only hope to enforce thresholds S=1+sqrt(eps) or larger');
end

if not(exist('maxSwaps','var')) || isempty(maxSwaps)
    maxSwaps=20*size(quad.X,2);
end

absdets = [];
maxnorms = [];
%optimization
swaps=0;
invcond=1;
indices = 1:length(quad.X); %helper for Matlab indexing
while(true)
    sym = symBasisFromQuadBasis(quad);
    absdets = [absdets abs(det(sym.X))];
    maxnorms = [maxnorms max(max(abs(sym.X)))];
    subplot(2,2,1);imagesc(abs(sym.X),[0 2*threshold]);
    subplot(2,2,2);plot(maxnorms);xlabel('iteration');ylabel('max |X|');
    subplot(2,2,3);semilogy(absdets);xlabel('iteration');ylabel('|det X|');
    drawnow;
    % TODO: use Natasa's heuristic instead of recomputing norms each time
    % column norms of C
    squaredColNorms = sum(abs(quad.X(quad.v,quad.v)).^2,1);
    [maxSquaredColNorm,maxi] = max(squaredColNorms);
    if maxSquaredColNorm > threshold
        % I can't find how to avoid this index juggling --federico
        indicesInV = indices(quad.v);
        maxi = indicesInV(maxi);
        [quad.X, quad.v, stepcond] = updateQuadBasisOut(quad.X, quad.v, maxi);
        swaps = swaps + 1;
    else
        % row norms of B
        squaredRowNorms = sum(abs(quad.X(~quad.v,~quad.v)).^2,2);
        [maxSquaredRowNorm,maxj] = max(squaredRowNorms);
        if maxSquaredRowNorm > threshold
            indicesNotInV = indices(~quad.v);
            maxj = indicesNotInV(maxj);
            [quad.X, quad.v, stepcond] = updateQuadBasisIn(quad.X, quad.v, maxj);
            swaps = swaps + 1;
        else
            [maxvec, maxis] = max(abs(quad.X(~quad.v,quad.v)));
            [maxval, maxj] = max(maxvec);
            if maxval > threshold
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
