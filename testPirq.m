function testPirq

reset(RandStream.getGlobalStream);

n=4;
m=2;

for tries=1:100
    U=randn(2*n+m,n+m);
    [Pi,R,Q,invcond]=pirq(U);
    
    assertVectorsAlmostEqual(abs(R*Q),abs(U(Pi,:)));
    assertVectorsAlmostEqual(Q*Q',eye(size(Q)));
    assertVectorsAlmostEqual(abs(Pi(1:n)-Pi(n+1:2*n))/n,ones(n,1));
end
