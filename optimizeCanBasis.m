function [can,invcond,swaps]=optimizeCanBasis(can,threshold,maxSwaps)
% given a canBasis, reduces it so that all elements are below a threshold
%
% [can,invcond,swaps]=optimizeCanBasis(can,threshold,maxswaps);
%
% for the meaning of invcond, see optimizeSymBasis
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

if isempty(can.X)
    invcond=1;swaps=0;
    return;
end

if not(exist('threshold','var')) || isempty(threshold)
    threshold=2;
elseif threshold<1+sqrt(eps(class(can.X)))
    error('cbrpack:thresholdTooSmall','you can only hope to enforce thresholds T=1+sqrt(eps) or larger');
end

if not(exist('maxSwaps','var')) || isempty(maxSwaps)
    maxSwaps=10*size(can.X,2);
end

assertEqual(length(can.p),size(can.X,2)+size(can.X,1));

invcond=1;
swaps=0;
while(swaps<maxSwaps)
    [maxvec, maxis]=max(abs(can.X));
    [maxval maxj]=max(maxvec);
    maxi=maxis(maxj);
    % the three lines above compute maxi,maxj=argmax(abs(X(i,j)))

    if maxval<threshold
        break;
    end

    [can.X,can.p,stepcond]=updateCanBasis(can.X,can.p,maxi,maxj);
    swaps=swaps+1;
    invcond=invcond*stepcond;
end
if swaps==maxSwaps
    warning('pgdoubling:stagnated','failed to produce a X with elements below the required threshold (obtained:%d, required:%d). Try running with a larger threshold.',maxval,threshold);
end
