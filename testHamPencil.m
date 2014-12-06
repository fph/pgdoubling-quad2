function testHamPencil
%tests pack/unpack symplectic pencils
reset(RandStream.getGlobalStream);

n=10;
for tries=1:100
    S=randomLagrangianSubspace(4*n);
    first=1:n;second=n+1:2*n;third=2*n+1:3*n;fourth=3*n+1:4*n;
    
    E = [S(first,:)' S(second,:)'];
    A = [-S(fourth,:)' S(third,:)'];
    
    assertElementsAlmostEqual(E*jay(2*n)*A',-A*jay(2*n)*E');
    
    sym=symBasisFromHamiltonianPencil(A,E);
    [A2,E2]=hamiltonianPencilFromSymBasis(sym);
    
    assertElementsAlmostEqual(subspace([A';E'],[A2';E2']),0);
end
