%reset(RandStream.getGlobalStream);
n = 10;
clear quad;
quad.v = [true(1,n) false(1,n)];
V = [randn(n,n-1) 1e6*randn(n,1)]; quad.X(quad.v,quad.v) = V;
W = [randn(n,n-1) 1e6*randn(n,1)]; quad.X(~quad.v,~quad.v) = W;
A = 1e6*randn(n); quad.X(~quad.v,quad.v) = A;
U = symplecticSubspaceFromQuadBasis(quad);

condU = cond(U)

quadOpt = optimizeQuadBasis(quad);
Uopt = symplecticSubspaceFromQuadBasis(quadOpt);
quad.v
quadOpt.v
distUopt = subspace(U,Uopt)

sym = symBasisFromQuadBasis(quad);
symOpt = optimizeSymBasis(sym);
UsymOpt = symplecticSubspaceFromSymBasis(symOpt);
distUsymOpt = subspace(U,UsymOpt)

X = symmetricMatrixFromSymBasis(sym);
Xopt = symmetricMatrixFromSymBasis(symBasisFromQuadBasis(quadOpt));
XsymOpt = symmetricMatrixFromSymBasis(symOpt);

%cholX = chol(X);
%cholXopt = chol(Xopt);
%cholXsymOpt = chol(XsymOpt);

errXopt = norm(X-Xopt)/norm(X)
errXsymOpt = norm(X-XsymOpt)/norm(X)

%norm(cholX-cholXopt)/norm(cholX)
%norm(cholX-cholXsymOpt)/norm(cholX)
