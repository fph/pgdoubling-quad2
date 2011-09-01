function U=symBasis2SymplecticSubspace(X,v)
% unpacks the symBasis representation of a subspace
%
% U=symBasis2SymplecticSubspace(X,v)
%
% given X,v representing a symplectic subspace, return a matrix spanning it.
%

U=[eye(size(X));X];
U=symplecticSwap(U,v,'T');
