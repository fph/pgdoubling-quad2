function testProductWithSparseL

reset(RandStream.getGlobalStream);

m=5;n=8;


for tries=1:100
    U=randn(m+n,n);
    
    can=canBasisFromSubspace(U);
    L=[eye(n) zeros(n,m); -can.X eye(m)];
    L(:,can.p)=L;
    tmp=L*U;
    assertVectorsAlmostEqual(tmp(n+1:end,:),zeros(m,n)); %checks that the constructed L works
    
    %multiplies with a random vector
    v=randn(m+n,2);
    assertVectorsAlmostEqual(productWithSparseL(can,v,'N'),L*v);
    assertVectorsAlmostEqual(productWithSparseL(can,v','T'),v'*L);
    assertVectorsAlmostEqual(productWithSparseL(can,v,'I'),L\v);
    assertVectorsAlmostEqual(productWithSparseL(can,v','IT'),v'/L);
end
