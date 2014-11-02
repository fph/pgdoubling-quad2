function testCanBasisFromSubspace

U=[1 1;1 1; 0 4];

can = canBasisFromSubspace(U);
U2 = subspaceFromCanBasis(can);

assertVectorsAlmostEqual(subspace(U,U2),0);

U=[1 1; 1 1; 1 1+1e-10];

assertWarningThrown(@() canBasisFromSubspace(U),'PGDoubling:illConditionedMatrix');

U=[4 0; 0 1; 10 20; 30 40];
[can,invcond,swaps]=canBasisFromSubspace(U,'permutation',[2 1 3 4]);
assertEqual(can.p,[2 1 3 4]);
assertEqual(can.X,[20 2.5;40 7.5]);
assertEqual(swaps,0);
assertEqual(can.origin,'subspace');

reset(RandStream.getGlobalStream);
for k=1:40
   U = randn(40,13).*exp(randn(40,13));
   can = canBasisFromSubspace(U, 'threshold', 1.5);
   assert(all(all(abs(can.X)<1.5)));
   newcan = canBasisFromSubspace(U);
   assertEqual(newcan.p,can.p);
   assert(all(all(abs(newcan.X)<1.5)));
end
