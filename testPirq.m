function testPirq

reset(RandStream.getGlobalStream);

n=4;

for tries=1:100
    U=randn(2*n,n);
    [v,R,Q,invcond]=pirq(U);
    
    assertVectorsAlmostEqual(R*Q,rowSwap(U,v,'N'));
    assertVectorsAlmostEqual(Q*Q',eye(size(Q)));
    assert(invcond>1e-8);
end
