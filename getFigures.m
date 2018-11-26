% getFigures.m

% Generate a random example of a quadBasis representation, run the
% optimization algorithm and produce the following figures:
%   * snapshots of |X| for iterations 1, 10, 20, final one (rng(0) selected
%   for the "stripes looks pretty" reason;
%   * max|x_{ij}| vs. iterations and |det X| vs. iterations 

rng(0);
n=30;
quad.X=randn(n);
quad.v=logical(randi(2,n,1)-1);
quad.v';
t = sum(quad.v);
fprintf(1,'C is %d x %d, A is %d x %d and B is %d x %d\n',t,t,n-t,t,n-t,n-t)

optimizeQuadBasis_figures(quad)
