function [A,E]=canBasis2Pencil(X,p)
% "unpack" a canonical basis returning a pencil
%
% [A,E]=canBasis2Subspace(X,p)
%
% inverse of pencil2CanBasis
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

EAt=canBasis2Subspace(X,p);

n=size(EAt,1)/2;

E=EAt(1:n,:)';
A=EAt(n+1:end,:)';