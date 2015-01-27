function testOptimizeQuadBasis

reset(RandStream.getGlobalStream);

for dt = [2 10 20]
    for odt = dt+[4 20 40]
        for trie = 1:30
            n=9;
            quad.X = randn(n) .* exp(4*randn(n)); quad.X = quad.X + quad.X';
            quad.v = randi(2,n,1)==1;
            quadOpt = optimizeQuadBasis(quad, dt, odt);
            U = symplecticSubspaceFromQuadBasis(quad);
            U2 = symplecticSubspaceFromQuadBasis(quadOpt);
            assertElementsAlmostEqual(subspace(U,U2),0);
            sym = symBasisFromQuadBasis(quadOpt);
            assertEqual(diag(sym.X)<dt,true(n,1));
            assertEqual(sym.X-diag(diag(sym.X))<odt,true(n,n));
        end
    end
end