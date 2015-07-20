reset(RandStream.getGlobalStream);
n = 7;
quad.v = [true(1,n) false(1,n)];
V = [randn(n,n-1) 1e12*randn(n,1)]; quad.X(quad.v,quad.v) = V;
W = [randn(n,n-1) 1e12*randn(n,1)]; quad.X(~quad.v,~quad.v) = W;
A = zeros(n); quad.X(~quad.v,quad.v) = A;
U = symplecticSubspaceFromQuadBasis(quad);

quadOpt = optimizeQuadBasis(quad);
Uopt = symplecticSubspaceFromQuadBasis(quadOpt);
quad.v
quadOpt.v
subspace(U,Uopt)

sym = symBasisFromQuadBasis(quad);
symOpt = optimizeSymBasis(sym);
UsymOpt = symplecticSubspaceFromSymBasis(symOpt);
subspace(U,UsymOpt)
