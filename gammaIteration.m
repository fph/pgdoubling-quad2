function gamma=gammaIteration(A1,B1,B2,C1,C2,D11,D12,D21,D22,tol,varargin);
% gamma iteration for H-infinity control
%
% gamma=gammaIteration(A1,B1,B2,C1,C2,D11,D12,D21,D22,tol,opts);
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

if size(varargin)==0
    o=Options('safer',true,'type','sign','maxSteps',100);
else
    o=Options(varargin{:});
end

[p1 m2]=size(D12);
[p2 m1]=size(D21);


eigsRH=eig([D11';D12']*[D11 D12], blkdiag(eye(m1),zeros(m2)));
eigsRJ=eig([D11;D21]*[D11' D21'], blkdiag(eye(p1),zeros(p2)));

gammaHatH=sqrt(max(eigsRH(isfinite(eigsRH))));
gammaHatJ=sqrt(max(eigsRJ(isfinite(eigsRJ))));

lowerBound=max(gammaHatH,gammaHatJ);
upperBound=inf;

gamma=norm([A1 B1 B2; C1 D11 D12; C2 D21 D22]); %half-assed starting value
reason='';

S=warning('off','cbrpack:illConditionedSubspace');
onCleanup(@() warning(S));

upperF=nan; %we want to use the secant method, therefore we need to keep track of the "y"s in the points upperBound and lowerBound
lowerF=nan;

iterations=0;

while(upperBound-lowerBound>tol)
    iterations=iterations+1;
    gammaOld=gamma;
    %gets the next value of gamma to try
    if upperBound==inf
        gamma=2*gamma;
    elseif not(isnan(lowerF+upperF)) && mod(iterations,5)==0
        %secant works well only in the very end, when we have accurate
        %eigenvalues, so we make one step every 5 with it
        gamma=lowerBound-lowerF*(upperBound-lowerBound)/(upperF-lowerF);
        assert(gamma>=lowerBound && gamma<=upperBound);
    else
        gamma=(upperBound+lowerBound)/2;
    end
    if gammaOld==gamma %we cannot get any more accurate than this apparently, ran into floating point max precision
        break;
    end
    
    if isempty(reason)
        %we are at the first iteration
        %TODO: print some pretty header
    else
        %completing line printed at the previous iteration
        fprintf('result: %s\n',reason);
    end
    
    fprintf('It %3d yLower=%e, yUpper=%e, xLower: %e, xUpper:%e, delta=%e, gamma tested=%e ',iterations,lowerF,upperF,lowerBound,upperBound,upperBound-lowerBound,gamma);
    
    [A,B,Q,R,S]=hInfinityControlPencil('J',gamma,A1,B1,B2,C1,C2,D11,D12,D21);
    [XJ YJ UJ VJ]=solveECARE(A,B,Q,R,S,o);
    k=checkECAREInvariantSubspaceResidual(A,B,Q,R,S,UJ);
    if not(k.pencilBackwardError<1e-10)
        lowerBound=gamma;
        lowerF=nan;
        reason='Could not solve the Riccati equation for X_J';
        continue;
    end
    
    [A,B,Q,R,S]=hInfinityControlPencil('H',gamma,A1,B1,B2,C1,C2,D11,D12,D21);
    [XH YH UH VH]=solveECARE(A,B,Q,R,S,o);
    k=checkECAREInvariantSubspaceResidual(A,B,Q,R,S,UH);
    if not(k.pencilBackwardError<1e-10)
        lowerBound=gamma;
        lowerF=nan;
        reason='Could not solve the Riccati equation for X_H';
        continue;
    end
    
%    %once we have found a stabilizable case, we stop updating expectedZeros
%    % since (1) the rank shouldn't change (2) computing the rank when
%    % close to the critical point is more ill-conditioned
    if(upperBound==inf) 
        nj=length(UJ)/2;
        nh=length(UH)/2;
        kj=rank(UJ(nj+1:end,:)'*UJ(1:nj,:),tol); %tol should not be too small
        kh=rank(UH(nh+1:end,:)'*UH(1:nh,:),tol);
        expectedZeros=nj+nh-kj-kh;
    end
    
    eigY=hInfinityEigenvaluesAlt(UH,UJ,gamma);
    eigY=real(eigY);
    [absoluteValues permutation]=sort(abs(eigY));
    eigY=eigY(permutation);
    assertElementsAlmostEqual(eigY(1:expectedZeros),zeros(expectedZeros,1));
    fValue=eigY(expectedZeros+1);
    if all(eigY(expectedZeros+1:end)>0)
        upperBound=gamma;
        upperF=fValue;
        reason='Stabilizable';
    else
        lowerBound=gamma;
        negativeIndex=find(eigY(expectedZeros+1:end)<0,1,'first');
        lowerF=eigY(expectedZeros+negativeIndex);
        reason='Not stabilizable';
    end
%    [expectedZeros eigY']
end
