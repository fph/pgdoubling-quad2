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

%S=warning('off','cbrpack:illConditionedMatrix');
%onCleanup(@() warning(S));
%T=warning('off','cbrpack:illConditionedMatrix');
%onCleanup(@() warning(T));

while(upperBound-lowerBound>tol)
    %gets the next value of gamma to try
    if upperBound==inf
        gamma=2*gamma;
    else
        %TODO: secant method!
        gamma=(upperBound+lowerBound)/2;
    end
    
    if isempty(reason)
        %we are at the first iteration
    else
        %completing line printed at the previous iteration
        fprintf('result: %s\n',reason);
    end
        
    fprintf('lower bound: %e, upper bound:%e, gamma tested=%e ',lowerBound,upperBound,gamma);
    
    [A,B,Q,R,S]=hInfinityControlPencil('J',gamma,A1,B1,B2,C1,C2,D11,D12,D21);
    [XJ YJ UJ VJ]=solveEcareSign(A,B,Q,R,S,'safer',true,'maxSteps',100);
    k=checkECAREInvariantSubspaceResidual(A,B,Q,R,S,UJ);
    if not(k.isGood)
        lowerBound=gamma;
        reason='Could not solve the Riccati equation for X_J';
    end
    
    [A,B,Q,R,S]=hInfinityControlPencil('H',gamma,A1,B1,B2,C1,C2,D11,D12,D21);
    [XH YH UH VH]=solveEcareSign(A,B,Q,R,S,'safer',true,'maxSteps',100);
    k=checkECAREInvariantSubspaceResidual(A,B,Q,R,S,UH);
    if not(k.isGood)
        lowerBound=gamma;
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
    if all(eigY(expectedZeros+1:end)>0)  
        upperBound=gamma;
        reason='Stabilizable';
    else
        lowerBound=gamma;
        reason='Not stabilizable';
    end
end
