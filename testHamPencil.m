function testHamPencil
%tests pack/unpack Hamiltonian pencils
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

%makes sure that the transformation choice is one such that v=all zeros
%corresponds to a positive-semidefinite sym.X when the signs are as in the
%standard control problem
[A, G, H] = carex(1); %not all CAREX problems work here; some have indefinite H, some do not have a symBasis with v=0
Ham = hamiltonian(A, G, H);
sym = symBasisFromHamiltonianPencil(Ham);
sym = symBasisFromSymBasis(sym,false(length(Ham),1));
assert(all(eig(sym.X)>-sqrt(eps)));