function testProductWithSparseL

reset(RandStream.getGlobalStream);

m=5;n=8;


for tries=1:100
    U=randn(m+n,n);
    
    [X,p]=subspace2CanBasis(U);
    L=[eye(n) zeros(n,m); -X eye(m)];
    L(:,p)=L;
    tmp=L*U;
    assertVectorsAlmostEqual(tmp(n+1:end,:),zeros(m,n)); %checks that the constructed L works
    
    %multiplies with a random vector
    v=randn(m+n,2);
    assertVectorsAlmostEqual(productWithSparseL(X,p,v,'N'),L*v);
    assertVectorsAlmostEqual(productWithSparseL(X,p,v','T'),v'*L);
    assertVectorsAlmostEqual(productWithSparseL(X,p,v,'I'),L\v);
    assertVectorsAlmostEqual(productWithSparseL(X,p,v','IT'),v'/L);
end
