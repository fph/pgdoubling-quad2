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
    w=randn(2,m+n);
    [Lv wL]=productWithSparseL(X,p,v,w);
    assertVectorsAlmostEqual(Lv,L*v);
    assertVectorsAlmostEqual(wL,w*L);
end
