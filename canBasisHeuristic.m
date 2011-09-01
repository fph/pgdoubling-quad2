function p=canBasisHeuristic(U,varargin)
% find a good starting permutation p for the canonical basis
% namely, U(p(1:size(U,2)),:) shouldn't be too ill-conditioned
%
% p=canBasisHeuristic(U)
%

[Q,R,p]=qr(U',0);
d=matgic.diagm(R);
if min(d)/max(d)<sqrt(eps(class(R)))
    warning('cbrpack:illConditionedSubspace','the provided subspace basis is ill-conditioned --- don''t expect much from the results');
end
