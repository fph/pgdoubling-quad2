function [X,p swaps invcond]=optimizeCanBasis(X,p,threshold,maxSwaps)
% given a canBasis, reduces it so that all elements are below a threshold
%
% [X,p swaps invcond]=optimizeCanBasis(X,p,threshold,maxswaps);
%
% for the meaning of invcond, see optimizeSymBasis

if not(exist('threshold','var')) || isempty(threshold)
    threshold=2;
elseif threshold<1+sqrt(eps(class(X)))
    error('cbrpack:thresholdTooSmall','you can only hope to enforce thresholds T=1+sqrt(eps) or larger');
end

if not(exist('maxSwaps','var')) || isempty(maxSwaps)
    maxSwaps=10*size(X,2);
end

assertEqual(length(p),size(X,2)+size(X,1));

invcond=1;
swaps=0;
while(swaps<maxSwaps)
    [maxvec, maxis]=max(abs(X));
    [maxval maxj]=max(maxvec);
    maxi=maxis(maxj);
    % the three lines above compute maxi,maxj=argmax(abs(X(i,j)))

    if maxval<threshold
        break;
    end

    [X,p,stepcond]=updateCanBasis(X,p,maxi,maxj);
    swaps=swaps+1;
    invcond=invcond*stepcond;
end
if swaps==maxSwaps
    warning('cbrpack:stagnated','failed to produce a X with elements below the required threshold (obtained:%d, required:%d). Try running with a larger threshold.',maxval,threshold);
end
