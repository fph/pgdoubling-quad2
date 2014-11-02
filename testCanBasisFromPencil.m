function testCanBasisFromPencil

reset(RandStream.getGlobalStream);

m=40;
n=25;
for tries=1:100
    A=5*randn(m,n)+5*1i*randn(m,n);
    E=5*randn(m,n)+5*1i*randn(m,n);
    
    can = canBasisFromPencil(A,E);
    [A2,E2] = pencilFromCanBasis(can);
    assertElementsAlmostEqual(subspace([E';A'],[E2';A2']),0); %checks that they are the same subspace
end

m=25;
n=40;
for tries=1:100
    A=5*randn(m,n)+5*1i*randn(m,n);
    E=5*randn(m,n)+5*1i*randn(m,n);
    
    can = canBasisFromPencil(A,E);
    [A2,E2] = pencilFromCanBasis(can);
    assertElementsAlmostEqual(subspace([E';A'],[E2';A2']),0); %checks that they are the same subspace
end
