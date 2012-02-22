function testHamiltonianMatrices

reset(RandStream.getGlobalStream);
n=6;
J=jay(n);

for tries=1:100
    X=randn(n);X=X+X';
    v=logical(randi(2,n,1)-1);
    [A,E]=symBasis2HamiltonianPencil(X,v);
    %generates Hamiltonian matrices
    assertVectorsAlmostEqual(A*J*E'+E*J*A',zeros(n));
    %pencil->symbasis->pencil
    v=logical(randi(2,n,1)-1);
    [A,E]=symBasis2HamiltonianPencil(X,v);
    assertVectorsAlmostEqual(A*J*E'+E*J*A',zeros(n));    
    M=randn(n);
    assertVectorsAlmostEqual(M*A*J*E'*M'+M*E*J*A'*M',zeros(n));
    [X,v]=hamiltonianPencil2SymBasis(M*A,M*E,v);
    [X,v]=optimizeSymBasis(X,v);
    [A2,E2]=symBasis2HamiltonianPencil(X,v);
    assertElementsAlmostEqual(subspace([A';E'],[A2';E2']),0);
end
