function [X,targetp]=canBasis2CanBasis(X,p,targetp)
% updates a canbasis to a new specified "target" p
%
% [X,targetp]=canBasis2CanBasis(X,p,targetp)
%
% equivalent to subspace2CanBasis(canBasis2Subspace(X,p),targetp);
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

%crappy implementation for now...
[X, targetp]=subspace2CanBasis(canBasis2Subspace(X,p),'permutation',targetp);
