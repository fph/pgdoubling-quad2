function testSymBasisFromSymBasis

n=5;

for trie = 1:100
    sym.v=logical(randi(2,n,1)-1);
    targetv=logical(randi(2,n,1)-1);
    sym.X=randn(n);sym.X=sym.X+sym.X';

    sym2=symBasisFromSymBasis(sym,targetv);
    assertEqual(sym2.v,targetv);
    symAlt = symBasisFromSymplecticSubspace(symplecticSubspaceFromSymBasis(sym),'swap',targetv);
    assertVectorsAlmostEqual(sym2.X,symAlt.X);
end