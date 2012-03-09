function p=canBasisHeuristic(U,varargin)
% find a good starting permutation p for the canonical basis
% namely, U(p(1:size(U,2)),:) shouldn't be too ill-conditioned
%
% p=canBasisHeuristic(U)
%

[Q,R,p]=qr(U',0);
