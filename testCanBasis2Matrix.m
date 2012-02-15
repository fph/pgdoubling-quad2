function testCanBasis2Matrix

reset(RandStream.getGlobalStream);

n=50;
r=47;
tries=50;
for i=1:tries
   M=randn(n,r)*5;
   [X,v]=matrix2CanBasis(M);
   [X,v]=optimizeCanBasis(X,v);
   N=canBasis2Matrix(X,v);
   assertVectorsAlmostEqual(M,N);
end
