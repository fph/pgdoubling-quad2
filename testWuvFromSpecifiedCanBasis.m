function testSpecifiedCanBasis2Wuv

n=10;
m=2;
r=3;

reset(RandStream.getGlobalStream);

for trie=1:100

wuv.u=randn(n,r);
wuv.w=randn(n,m);
wuv.v=randn(n,r);

U=[eye(n)-wuv.u*wuv.v';wuv.w'];

can=canBasisFromWuv(wuv);
can2=specifiedCanBasisFromSubspace(U,can.p);
assertVectorsAlmostEqual(can.X,can2.X);
end
