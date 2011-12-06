function testSymPencil
%tests pack/unpack symplectic pencils
reset(RandStream.getDefaultStream);

n=10;
for tries=1:100
S=randomLagrangianSubspace(4*n);
first=1:n;second=n+1:2*n;third=2*n+1:3*n;fourth=3*n+1:4*n;

U=[S(first,:)' S(third,:)'];
L=[S(fourth,:)' S(second,:)'];

assertElementsAlmostEqual(U*matgic.jay(2*n)*U',L*matgic.jay(2*n)*L');

[X,v]=symplecticPencil2SymBasis(L,U);
[L2,U2]=symBasis2SymplecticPencil(X,v);

assertElementsAlmostEqual(subspace([L';U'],[L2';U2']),0);
end
