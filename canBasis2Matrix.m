function X=canBasis2Matrix(X,p);
% turns a canbasis into a matrix
% 
% M=canBasis2Matrix(X,p);
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling


% in fact, the matrix representation is a particular canBasis, so we just
% have to get it

n=size(X,2);
in=find(p(n+1:end)<=n);
out=find(p(1:n)>n);
[X,p]=updateCanBasis(X,p,in,out);
X(p(n+1:end)-n,p(1:n))=X;
X=X';
