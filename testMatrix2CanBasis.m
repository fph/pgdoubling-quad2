function testMatrix2CanBasis

reset(RandStream.getGlobalStream);

m=40;
n=25;
for tries=1:100
    M=5*randn(m,n);
    
    [X,p]=matrix2CanBasis(M);
    assertVectorsAlmostEqual(M,canBasis2Matrix(X,p));
end

m=25;
n=40;
for tries=1:100
    M=5*randn(m,n);
    
    [X,p]=matrix2CanBasis(M);
    assertVectorsAlmostEqual(M,canBasis2Matrix(X,p));
end
