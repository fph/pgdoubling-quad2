function [X,v,swaps]=optimizeSymBasis(X,v,diagonalThreshold,offdiagonalThreshold,maxSwaps)
% given a symBasis, reduces it so that all elements are below a threshold
%
% [X,v,swaps]=optimizeSymBasis(X,v,diagonalThreshold,offdiagonalThreshold,maxSwaps)
%


if not(exist('diagonalThreshold','var')) || isempty(diagonalThreshold)
    diagonalThreshold=2;
end
    
if not(exist('odddiagonalThreshold','var')) || isempty(offdiagonalThreshold)
    offdiagonalThreshold=sqrt(5+diagonalThreshold^2);
end

if diagonalThreshold<1+sqrt(eps(class(X)))
    error('cbrpack:thresholdTooSmall','you can only hope to enforce thresholds S=1+sqrt(eps) or larger');
end

if offdiagonalThreshold<sqrt(1+diagonalThreshold^2)+sqrt(eps(class(X)))
    error('cbrpack:thresholdTooSmall','you can only hope to enforce thresholds T=sqrt(1+S^2) or larger');
end

if not(exist('maxSwaps','var')) || isempty(maxSwaps)
    maxSwaps=10*size(X,2);
end

%optimization
swaps=0;
while(swaps<maxSwaps)
    [maxvec, maxis]=max(abs(X-diag(diag(X))));
    [maxval maxj]=max(maxvec);
    maxi=maxis(maxj);
    % the three lines above compute maxi,maxj=argmax(abs(X(i,j)))

    [maxdiag maxdiagPos]=max(diag(X));
    if maxdiag>diagonalThreshold
        [X,v]=updateSymBasis(X,v,maxdiagPos);
        swaps=swaps+1;
    elseif maxval>offdiagonalThreshold
        [X,v]=updateSymBasis(X,v,[maxi maxj]);
        swaps=swaps+2;
    else
        break;
    end    
end

if swaps==maxSwaps
    warning('cbrpack:stagnated','failed to produce a X with elements below the required threshold (obtained:%d, required:%d). Try running with a larger threshold.',maxval,threshold);
end