function testCanBasisFromCanBasis

m=3;
n=8;

reset(RandStream.getGlobalStream);

for tries=1:100
    can.p=randperm(n+m);
    can.X=randn(m,n);
    can.origin='subspace';
    targetp=randperm(n+m);
    
    newcan=canBasisFromCanBasis(can,targetp);
    assertEqual(newcan.p,targetp);
    can2 = canBasisFromSubspace(subspaceFromCanBasis(can),'permutation',targetp);
    assertVectorsAlmostEqual(can2.X, newcan.X);
end
