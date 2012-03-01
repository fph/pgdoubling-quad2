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

% start reducing the last block to triangular form, then invokes pirq
[Q,R]=qr(U(2*n+1:end,:)');
d=abs(diag(R));
if min(d)/max(d)<sqrt(eps(class(U)))
    warning('cbrpack:illConditionedThirdBlock','the extended part of the "extended Lagrangian" subspace is almost singular --- you are about to invert an ill-conditioned matrix');
end

U=U*Q; %now U3 should be block lower triangular
[v,R2,Q2,invcond]=pirq(U(1:2*n,m+1:end));

if invcond<sqrt(eps(class(U)))
    warning('cbrpack:illConditionedSubspace','the provided subspace basis is ill-conditioned --- don''t expect much from the results');
end
