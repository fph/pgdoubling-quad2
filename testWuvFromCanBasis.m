function testWuvFromCanBasis

reset(RandStream.getGlobalStream);

n=5;m=16;

for trie=1:100
   can.X=randn(n,m);
   can.p=randperm(n+m);
   can.origin='matrix';
   M=matrixFromCanBasis(can);
   wuv = wuvFromCanBasis(can);
   assertElementsAlmostEqual(M,wuv.w'/(eye(m)-wuv.u*wuv.v'));
   assertElementsAlmostEqual(eye(size(wuv.u,2))-wuv.v'*wuv.u,wuv.ImvTu);
end
