function U=canBasis2Subspace(X,p)
% "unpack" a canonical basis returning a matrix spanning the subspace
%
% U=canBasis2Subspace(X,p)
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling


[m n]=size(X);
p=p(:);

if not(length(p)==m+n)
    error('cbrpack:wrongPermutationLength','length of the permutation vector does not match size of X (length(p)=%d,expected %d)',length(p),sum(size(X)));
end

n=size(X,2);

U=zeros(m+n,n);
U(p(1:n),:)=eye(n);
U(p(n+1:end),:)=X;
