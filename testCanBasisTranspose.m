function testCanBasisTranspose

reset(RandStream.getGlobalStream);

n=12;m=7;

for trie=1:100
    X=randn(n,m);
    p=randperm(n+m);
    [Xt,pt]=canBasisTranspose(X,p);
    M=canBasis2Matrix(X,p);
    Mt=canBasis2Matrix(Xt,pt);
    assertVectorsAlmostEqual(M',Mt);
end

n=4;m=13;

for trie=1:100
    X=randn(n,m);
    p=randperm(n+m);
    [Xt,pt]=canBasisTranspose(X,p);
    M=canBasis2Matrix(X,p);
    Mt=canBasis2Matrix(Xt,pt);
    assertVectorsAlmostEqual(M',Mt);
end
