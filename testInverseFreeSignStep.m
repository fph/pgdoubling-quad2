function testInverseFreeSignStep
n=6;
X=randn(n);X=X+X';v=logical(randi(2,n,1)-1);

[A,E]=symBasis2HamiltonianPencil(X,v);
H=inv(E)*A;
Hnew=1/2*(H+inv(H));

[Xnew,vnew,w,swaps1,swaps2,res]=inverseFreeSignStep(X,v);

[Anew,Enew]=symBasis2HamiltonianPencil(Xnew,vnew);
assertElementsAlmostEqual(0,subspace([Anew';Enew'],[Hnew';eye(size(Hnew))]));
