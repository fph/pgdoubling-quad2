function [X,v]=defectCorrectionCARE(A,G,Q,X,v,steps)
%implements defect correction for a CARE
%
% [X,v]=defectCorrectionCARE(A,G,Q,X,v)

H=hamiltonian(A,G,Q);

%switches to Pi'*H*Pi, whose invariant subspace is [I;X]

H=rowSwap(rowSwap(H,v,'T')',v,'T')';

%assertElementsAlmostEqual(0,subspace(H*[eye(size(X));X],[eye(size(X));X]));

[A,G,Q]=Hamiltonian2RiccatiCoefficients(H);

AA=A-G*X;
QQ=Q+X*A+A'*X-X*G*X;

norm(AA),norm(QQ),norm(G)

DeltaX=solveCARE(AA,G,QQ,'verbose',true,'maxSteps',steps,'minSteps',steps);
norm(X),norm(DeltaX)
X=X+DeltaX;
