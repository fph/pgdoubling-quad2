function testWuv2SpecifiedCanBasis

n=10;
m=2;
r=3;

reset(RandStream.getGlobalStream);

for trie=1:100

u=randn(n,r);
w=randn(n,m);
v=randn(n,r);

U=[eye(n)-u*v';w'];

[X,p]=wuv2CanBasis(w,u,v);
X2=subspace2SpecifiedCanBasis(U,p);
assertVectorsAlmostEqual(X,X2);
end
