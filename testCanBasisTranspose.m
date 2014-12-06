function testCanBasisTranspose

reset(RandStream.getGlobalStream);

n=12;m=7;

for trie=1:100
    can.X=randn(n,m);
    can.p=randperm(n+m);
    can.origin='matrix';
    cant=canBasisTranspose(can);
    M=matrixFromCanBasis(can);
    Mt=matrixFromCanBasis(cant);
    assertVectorsAlmostEqual(M',Mt);
end

n=4;m=13;

for trie=1:100
    can.X=randn(n,m);
    can.p=randperm(n+m);
    can.origin='matrix';
    cant=canBasisTranspose(can);
    M=matrixFromCanBasis(can);
    Mt=matrixFromCanBasis(cant);
    assertVectorsAlmostEqual(M',Mt);
end
