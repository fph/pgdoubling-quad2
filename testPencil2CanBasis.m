function testCanBasis2Subspace


reset(RandStream.getGlobalStream);

m=40;
n=25;
for tries=1:100
    A=5*randn(m,n)+5*1i*randn(m,n);
    E=5*randn(m,n)+5*1i*randn(m,n);
    
    [X,p]=pencil2CanBasis(A,E);
    [A2,E2]=canBasis2Pencil(X,p);
    assertElementsAlmostEqual(subspace([E';A'],[E2';A2']),0); %checks that they are the same subspace
end

m=25;
n=40;
for tries=1:100
    A=5*randn(m,n)+5*1i*randn(m,n);
    E=5*randn(m,n)+5*1i*randn(m,n);
    
    [X,p]=pencil2CanBasis(A,E);
    [A2,E2]=canBasis2Pencil(X,p);
    assertElementsAlmostEqual(subspace([E';A'],[E2';A2']),0); %checks that they are the same subspace
end
