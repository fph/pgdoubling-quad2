function [v invcond X]=symBasisHeuristic(U)
% find a good starting row swap v for a symplectic basis
% namely, the top half of Pi_v^T(U) shouldn't be too ill-conditioned
%
% [v invcond X]=symBasisHeuristic(U)
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

[v,R,Q,invcond]=pirq(U);

n=length(U);

if nargout>2
    [X invcond]=rightLinSolve(R(n+1:end,:),R(1:n,:));
end

if invcond<sqrt(eps(class(R)))
    warning('cbrpack:illConditionedSubspace','the provided subspace basis is ill-conditioned --- don''t expect much from the results');
end
