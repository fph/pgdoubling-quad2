function testHamiltonianMatrices

reset(RandStream.getGlobalStream);
n=6;
J=jay(n);

for tries=1:100
    sym.X=randn(n);sym.X=sym.X+sym.X';
    sym.v=logical(randi(2,n,1)-1);
    sym.origin='hamiltonianPencil';
    [A,E]=hamiltonianPencilFromSymBasis(sym);
    %generates Hamiltonian matrices
    assertVectorsAlmostEqual(A*J*E'+E*J*A',zeros(n));
    %pencil->symbasis->pencil
    sym.v=logical(randi(2,n,1)-1);
    [A,E]=hamiltonianPencilFromSymBasis(sym);
    assertVectorsAlmostEqual(A*J*E'+E*J*A',zeros(n));    
    M=randn(n);
    assertVectorsAlmostEqual(M*A*J*E'*M'+M*E*J*A'*M',zeros(n));
    sym=symBasisFromHamiltonianPencil(M*A,M*E,'initialSwap',sym.v);
    [A2,E2]=hamiltonianPencilFromSymBasis(sym);
    assertElementsAlmostEqual(subspace([A';E'],[A2';E2']),0);
end
