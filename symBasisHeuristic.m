function [v invcond]=symBasisHeuristic(U)
% find a good starting row swap v for a symplectic basis
% namely, the top half of Pi_v^T(U) shouldn't be too ill-conditioned
%
% [v invcond]=symBasisHeuristic(U)
%

[v,R,q,invcond]=pirq(U);

if invcond<sqrt(eps(class(R)))
    warning('cbrpack:illConditionedSubspace','the provided subspace basis is ill-conditioned --- don''t expect much from the results');
end
