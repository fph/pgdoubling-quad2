function testCanBasisFromMatrix

reset(RandStream.getGlobalStream);

m=40;
n=25;
for tries=1:100
    M=5*randn(m,n);
    
    can = canBasisFromMatrix(M);
    assertEqual(can.origin,'matrix');
    assertVectorsAlmostEqual(M,matrixFromCanBasis(can));
end

m=25;
n=40;
for tries=1:100
    M=5*randn(m,n);
    
    can = canBasisFromMatrix(M);
    assertVectorsAlmostEqual(M,matrixFromCanBasis(can));
end
