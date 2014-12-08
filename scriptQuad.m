n=7;
k=2;
reset(RandStream.getGlobalStream);
quad.X = randn(n);
quad.v = [true(k,1); false(n-k,1)];
quad.X(quad.v,~quad.v)=0;
% for now we think M = [-A'*A 0; B C*C'] and we store X=[A 0; B C]
% We'll set v(i)=true iff i belongs to the first block
% the lower triangle of X is used only.

oldsym.X = matrixFromQuadBasis(quad);
oldsym.v = quad.v;
targetv = quad.v; targetv(i) = ~targetv(i);
newsym = symBasisFromSymBasis(oldsym, targetv);

oldquad=quad;

i = k+1; %entry to PPT on
assert(quad.v(i)==false); %for now we assume we go in this direction

oldv = quad.v;
quad.v(i) = true;
[w, beta, c] = gallery('house', quad.X(i,~oldv)');
% apply Householder transformation
quad.X(~oldv,~oldv) = quad.X(~oldv,~oldv) - quad.X(~oldv,~oldv)*w*beta*w';
% update B block
% interesting observation: what follows is a PPT on the B block
quad.X(~quad.v,i) = quad.X(~quad.v,i) / c;
quad.X(~quad.v,oldv) = quad.X(~quad.v,oldv) - quad.X(~quad.v,i) * quad.X(i,oldv);
% update A block
quad.X(oldv,i) = 0;
quad.X(i, oldv) = - quad.X(i, oldv) / c;
quad.X(i,i) = 1/c;

norm(matrixFromQuadBasis(quad) - newsym.X)
