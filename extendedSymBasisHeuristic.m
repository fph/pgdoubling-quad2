function v=extendedSymBasisHeuristic(U)
% computes a "reasonable" value of v for an extended Lagrangian subspace
%
% v=symBasisHeuristic(U)
%
% uses QRP (just a crappy heuristic)
%
% U is 2n+m x n+m

[a b]=size(U);
n=a-b;
m=b-n;

[Q,R,E]=qr(U',0);
d=abs(diag(R));

if min(d)/max(d)<sqrt(eps(class(R)))
    warning('cbrpack:illConditionedSubspace','the provided subspace basis is ill-conditioned --- don''t expect much from the results');
end

invE=1:length(E);
invE(E)=invE;

v=(invE(1:n)>invE(n+1:2*n));
