function testSymBasis2symBasis

n=5;

v=logical(randi(2,n,1)-1);
targetv=logical(randi(2,n,1)-1);
X=randn(n);X=X+X';

[Xup,vup]=symBasis2SymBasis(X,v,targetv);
assertEqual(vup,targetv);

assertVectorsAlmostEqual(Xup,symplecticSubspace2SymBasis(symBasis2SymplecticSubspace(X,v),'swap',targetv));
