function testCanBasis2wuv

reset(RandStream.getGlobalStream);

n=5;m=16;

for trie=1:100
   X=randn(n,m);
   p=randperm(n+m);
   M=canBasis2Matrix(X,p);
   [w u v ImvTu]=canBasis2Wuv(X,p);
   assertElementsAlmostEqual(M,w'/(eye(m)-u*v'));
   assertElementsAlmostEqual(eye(size(u,2))-v'*u,ImvTu);
end