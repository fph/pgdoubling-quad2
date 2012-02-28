function [v invcond]=extendedSymBasisHeuristicPaper(U)
% computes a "reasonable" value of v for an extended Lagrangian subspace
%
% [v invcond]=symBasisHeuristic(U)
%
% computes a "restricted QRP" 
%
% U is 2n+m x n+m, 1:n can be swapped with n+1:2*n

[a b]=size(U);
n=a-b;
m=b-n;

% computes a "restricted" (QRP)' factorization

[Pi R Q invcond]=pirq(U);

invPi=1:length(Pi);
invPi(Pi)=invPi;

v=(invPi(1:n)>invPi(n+1:2*n));

if invcond<sqrt(eps(class(R)))
    warning('cbrpack:illConditionedSubspace','the provided subspace basis is ill-conditioned --- don''t expect much from the results');
end
