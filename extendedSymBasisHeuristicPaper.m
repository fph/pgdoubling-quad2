function [v invcond]=extendedSymBasisHeuristicPaper(U)
% computes a "reasonable" value of v for an extended Lagrangian subspace
%
% [v invcond]=symBasisHeuristic(U)
%
% computes a "restricted QRP" 
%
% U is 2n+m x n+m, 1:n can be swapped with n+1:2*n
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

[a b]=size(U);
n=a-b;
m=b-n;

% start reducing the last block to triangular form, then invokes pirq
[Q,R]=qr(U(2*n+1:end,:)');
if m>1 %the badly-designed Matlab diag does not do what I mean for m=1 :(
    d=abs(diag(R));
else
    d=abs(R(1,1));
end
if min(d)/max(d)<sqrt(eps(class(U)))
    warning('cbrpack:illConditionedThirdBlock','the extended part of the "extended Lagrangian" subspace is almost singular --- you are about to invert an ill-conditioned matrix');
end

U=U*Q; %now U3 should be block lower triangular
[v,R2,Q2,invcond]=pirq(U(1:2*n,m+1:end));

if invcond<sqrt(eps(class(U)))
    warning('cbrpack:illConditionedSubspace','the provided subspace basis is ill-conditioned --- don''t expect much from the results');
end
