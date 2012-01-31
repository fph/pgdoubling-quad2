function testDoublingStep

n=10;
X=randn(n);X=(X+X')/2;
v=logical(randi(2,n,1)-1);
[L,U]=symBasis2SymplecticPencil(X,v);
[X2,v2]=doublingStep(X,v,[]);
M=U\L;
[XX,vv]=symplecticPencil2SymBasis(M^2,eye(n));

assertElementsAlmostEqual(subspace(symBasis2SymplecticSubspace(X2,v2),symBasis2SymplecticSubspace(XX,vv)),0);
