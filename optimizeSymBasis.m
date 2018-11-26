function [symb,invcond,swaps]=optimizeSymBasis(symb,diagonalThreshold,offDiagonalThreshold,maxSwaps)
% given a symBasis, reduces it so that all elements are below a threshold
%
% [symb,invcond,swaps]=optimizeSymBasis(symb,diagonalThreshold,offdiagonalThreshold,maxSwaps)
%
% invcond is the inverse condition number of all the performed
% transformations composed. If it is larger than some moderate value,
% consider recomputing the basis
%
% TODO: the condition number probably isn't the right thing --- if X=[1e12
% 1; 1 1], then invcond is huge but the computation is perfectly
% conditioned
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

if not(exist('diagonalThreshold','var')) || isempty(diagonalThreshold)
    diagonalThreshold=2;
end

if not(exist('offDiagonalThreshold','var')) || isempty(offDiagonalThreshold)
    offDiagonalThreshold=sqrt(5+diagonalThreshold^2);
end

if diagonalThreshold<1+sqrt(eps(class(symb.X)))
    error('PGDoubling:thresholdTooSmall','you can only hope to enforce thresholds S=1+sqrt(eps) or larger');
end

if offDiagonalThreshold<sqrt(1+diagonalThreshold^2)+sqrt(eps(class(symb.X)))
    error('PGDoubling:thresholdTooSmall','you can only hope to enforce thresholds T=sqrt(1+S^2) or larger');
end

if not(exist('maxSwaps','var')) || isempty(maxSwaps)
    maxSwaps=20*size(symb.X,2);
end

%optimization
swaps=0;
invcond=1;
while(true)
    [maxvec, maxis] = max(abs(symb.X-diag(diag(symb.X))));
    [maxval, maxj] = max(maxvec);
    maxi = maxis(maxj);
    % the three lines above compute maxi,maxj = argmax(abs(X(i,j)))
    
    [maxdiag, maxdiagPos] = max(abs(diag(symb.X)));
    if maxdiag > diagonalThreshold
        [symb.X,symb.v,stepcond] = updateSymBasis(symb.X,symb.v,maxdiagPos);
        swaps = swaps+1;
    elseif maxval > offDiagonalThreshold
        [symb.X,symb.v,stepcond] = updateSymBasis(symb.X,symb.v,[maxi maxj]);
        swaps = swaps+2;
    else
        break;
    end
    invcond = invcond*stepcond;
    if swaps >= maxSwaps
        warning('PGDoubling:stagnated','failed to produce a X with elements below the required threshold. Try running with a larger threshold.');
        break
    end
end
