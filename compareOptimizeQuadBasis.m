function testOptimizeQuadBasis

reset(RandStream.getGlobalStream);

for threshold = [1.1 2 10 20]
    for trie = 1:25
        n=100;
%        quad.X = randn(n) .* exp(4*randn(n));
        quad.X = randn(n) .* exp(3*randn(n));
        quad.v = randi(2,n,1)==1;
        [quadOpt,invcond2] = optimizeQuadBasis(quad, threshold);
        U = symplecticSubspaceFromQuadBasis(quad);
        U2 = symplecticSubspaceFromQuadBasis(quadOpt);
        
        %%% temporary code to check conditioning of two different
        %%% approaches
        sym = symBasisFromQuadBasis(quad);
        [symOpt,invcond3] = optimizeSymBasis(sym);
        U3 = symplecticSubspaceFromSymBasis(symOpt);
        quadRecomputed = quadBasisFromQuadBasis(quad,quadOpt.v);
        U4 = symplecticSubspaceFromQuadBasis(quadRecomputed);
        symRecomputed = symBasisFromSymBasis(sym,symOpt.v);
        U5 = symplecticSubspaceFromSymBasis(symRecomputed);
        [subspace(U,U2),subspace(U,U3),subspace(U,U4),subspace(U,U5),invcond2,invcond3,cond(U)]

        %%%
        
        assertElementsAlmostEqual(subspace(U,U2),0);
        sym = symBasisFromQuadBasis(quadOpt);
        assertEqual(sym.X<threshold,true(size(sym.X)));
    end
end