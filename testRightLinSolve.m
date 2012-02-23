function testRightLinSolve

reset(RandStream.getGlobalStream);
n=5;
m=2;

b=randn(m,n);
A=randn(n);

[X,c]=rightLinSolve(b,A);

assertVectorsAlmostEqual(X*A,b);
