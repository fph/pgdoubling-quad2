function [X,targetp]=canBasis2canBasis(X,p,targetp)
% updates a canbasis to a new specified "target" p
%
% [X,targetp]=canBasis2canBasis(X,p,targetp)
%
% equivalent to subspace2CanBasis(canBasis2Subspace(X,p),targetp);

%crappy implementation for now...
[X, targetp]=subspace2CanBasis(canBasis2Subspace(X,p),'permutation',targetp);
