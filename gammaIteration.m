function gamma=gammaIteration(A1,B1,B2,C1,C2,D11,D12,D21,D22,tol);
% gamma iteration for H-infinity control
%
% gamma=gammaIteration(A1,B1,B2,C1,C2,D11,D12,D21,D22,tol);

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
    %gets the next value of gamma to try
    if upperBound==inf
        gamma=2*gamma;
    elseif not(isnan(lowerF+upperF)) && mod(iterations,4)~=0
        %the secant method can still be dramatically slower than bisection,
        %so we make sure we try a little bit of both and force a bisection
        %step every 4 iterations
        gamma=lowerBound-lowerF*(upperBound-lowerBound)/(upperF-lowerF);
        assert(gamma>lowerBound && gamma<upperBound);
    else
        gamma=(upperBound+lowerBound)/2;
    end
    
    if isempty(reason)
        %we are at the first iteration
        %TODO: print some pretty header
    else
        %completing line printed at the previous iteration
        fprintf('result: %s\n',reason);
    end
        
    fprintf('It %3d yLower=%e, yUpper=%e, xLower: %e, xUpper:%e, gamma tested=%e ',iterations,lowerF,upperF,lowerBound,upperBound,gamma);
    
    [A,B,Q,R,S]=hInfinityControlPencil('J',gamma,A1,B1,B2,C1,C2,D11,D12,D21);
    [XJ YJ UJ VJ]=solveEcareSign(A,B,Q,R,S,'safer',true,'maxSteps',100);
    k=checkECAREInvariantSubspaceResidual(A,B,Q,R,S,UJ);
    if not(k.isGood)
        lowerBound=gamma;
        lowerF=nan;
        reason='Could not solve the Riccati equation for X_J';
    end
    
    [A,B,Q,R,S]=hInfinityControlPencil('H',gamma,A1,B1,B2,C1,C2,D11,D12,D21);
    [XH YH UH VH]=solveEcareSign(A,B,Q,R,S,'safer',true,'maxSteps',100);
    k=checkECAREInvariantSubspaceResidual(A,B,Q,R,S,UH);
    if not(k.isGood)
        lowerBound=gamma;
        lowerF=nan;
        reason='Could not solve the Riccati equation for X_H';
    end

    nj=length(UJ)/2;
    nh=length(UH)/2;
    kj=rank(UJ(nj+1:end,:)'*UJ(1:nj,:),tol);
    kh=rank(UH(nh+1:end,:)'*UH(1:nh,:),tol);

    eigY=hInfinityEigenvaluesAlt(UH,UJ,gamma);
    [absoluteValues permutation]=sort(abs(eigY));
    eigY=eigY(permutation);
    expectedZeros=nj+nh-kj-kh;
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
end
