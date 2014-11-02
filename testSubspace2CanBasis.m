function testSubspace2CanBasis

U=[1 1;1 1; 0 4];

[X,p]=subspace2CanBasis(U);
U2=canBasis2Subspace(X,p);

assertVectorsAlmostEqual(subspace(U,U2),0);

U=[1 1; 1 1; 1 1+1e-10];

assertWarningThrown(@() subspace2CanBasis(U),'cbrpack:badSubspace');

U=[4 0; 0 1; 10 20; 30 40];
[X,p]=subspace2CanBasis(U,'permutation',[2 1 3 4]);
assertEqual(p,[2 1 3 4]);
assertEqual(X,[20 2.5;40 7.5]);

reset(RandStream.getGlobalStream);
for k=1:40
   U=randn(40,13).*exp(randn(40,13));
   [X,p]=subspace2CanBasis(U,'threshold',1.5);
   assert(all(all(abs(X)<1.5)));
   [X,pnew]=subspace2CanBasis(U);
   assertEqual(pnew,p);
   assert(all(all(abs(X)<1.5)));
end
