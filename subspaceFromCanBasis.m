function U=subspaceFromCanBasis(can)
% "unpack" a canonical basis returning a matrix spanning the subspace
%
% U = canBasisFromSubspace(can)
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

if not(strcmp(can.origin,'subspace'))
    warning('PGDoubling:wrongOrigin','You are converting to a subspace a canbasis that did not originate from a subspace');
end

[m, n] = size(can.X);
p = can.p(:);

if not(length(p) == m+n)
    error('PGDoubling:wrongPermutationLength','the length of the permutation vector does not match the size of X (length(p)=%d,expected %d)',length(can.p),sum(size(can.X)));
end

n = size(can.X,2);

U = zeros(m+n,n);
U(p(1:n),:) = eye(n);
U(p(n+1:end),:) = can.X;
