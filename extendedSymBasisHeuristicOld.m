function [v invcond]=extendedSymBasisHeuristic(U)
% computes a "reasonable" value of v for an extended Lagrangian subspace
%
% [v invcond]=symBasisHeuristic(U)
%
% uses QRP (just a crappy heuristic)
%
% U is 2n+m x n+m
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

[a b]=size(U);
n=a-b;
m=b-n;

[Q,R,E]=qr(U',0);
d=abs(diag(R));

invcond=min(d)/max(d);

if invcond<sqrt(eps(class(R)))
    warning('cbrpack:illConditionedSubspace','the provided subspace basis is ill-conditioned --- don''t expect much from the results');
end

invE=1:length(E);
invE(E)=invE;

v=(invE(1:n)>invE(n+1:2*n));
