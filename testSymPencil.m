function testSymPencil
%tests pack/unpack symplectic pencils
reset(RandStream.getGlobalStream);

n=10;
for tries=1:100
S=randomLagrangianSubspace(4*n);
first=1:n;second=n+1:2*n;third=2*n+1:3*n;fourth=3*n+1:4*n;

U=[S(first,:)' S(third,:)'];
L=[S(fourth,:)' S(second,:)'];

assertElementsAlmostEqual(U*jay(2*n)*U',L*jay(2*n)*L');

sym=symBasisFromSymplecticPencil(L,U);
[L2,U2]=symplecticPencilFromSymBasis(sym);

assertElementsAlmostEqual(subspace([L';U'],[L2';U2']),0);
end
