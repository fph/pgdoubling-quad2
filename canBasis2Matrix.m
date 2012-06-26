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

[m n]=size(X);
X=canBasis2CanBasis(X,p,1:m+n);
