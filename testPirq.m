function testPirq

reset(RandStream.getGlobalStream);

n=4;
m=2;

for tries=1:100
    U=randn(2*n+m,n+m);
    [Pi,R,Q]=pirq(U);
    
    assertVectorsAlmostEqual(R*Q,U(Pi,:));
    assertVectorsAlmostEqual(Q*Q',eye(size(Q)));
    assertVectorsAlmostEqual(mod(Pi(1:n),n), mod(Pi(n+1:2*n),n));
end
