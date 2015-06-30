function testOptimizeQuadBasis

reset(RandStream.getGlobalStream);

for threshold = [1.1 2 10 20]
    for trie = 1:25
        n=15;
        quad.X = randn(n) .* exp(4*randn(n));
        quad.v = randi(2,n,1)==1;
        quadOpt = optimizeQuadBasis(quad, threshold);
        U = symplecticSubspaceFromQuadBasis(quad);
        U2 = symplecticSubspaceFromQuadBasis(quadOpt);
        assertElementsAlmostEqual(subspace(U,U2),0);
        sym = symBasisFromQuadBasis(quadOpt);
        assertEqual(sym.X<threshold,true(size(sym.X)));
    end
end