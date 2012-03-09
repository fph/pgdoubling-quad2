function testCanBasis2CanBasis

m=3;
n=8;

reset(RandStream.getGlobalStream);

for tries=1:100
    p=randperm(n+m);
    targetp=randperm(n+m);
    
    X=randn(m,n);
    
    [Xup,pup]=canBasis2canBasis(X,p,targetp);
    assertEqual(pup,targetp);
    
    assertVectorsAlmostEqual(Xup,subspace2CanBasis(canBasis2Subspace(X,p),'permutation',targetp));
end
