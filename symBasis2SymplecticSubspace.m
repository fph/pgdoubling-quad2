function U=symBasis2SymplecticSubspace(X,v)
% unpacks the symBasis representation of a subspace
%
% U=symBasis2SymplecticSubspace(X,v)
%
% given X,v representing a symplectic subspace, return a matrix spanning it.
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

U=[eye(size(X));X];
U=rowSwap(U,v,'T');
