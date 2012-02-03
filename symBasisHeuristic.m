function v=symBasisHeuristic(U,varargin)
% find a good starting row swap v for a symplectic basis
% namely, the top half of Pi_v^T(U) shouldn't be too ill-conditioned
%
% v=symBasisHeuristic(U)
%

n=size(U,2);
if(size(U,1)~=2*n)
    error('cbrpack:oddSize','the input matrix must have an even number of rows');
end

[Q,R,E]=qr(U',0);
d=abs(matgic.diagm(R));

if min(d)/max(d)<sqrt(eps(class(R)))
    warning('cbrpack:illConditionedSubspace','the provided subspace basis is ill-conditioned --- don''t expect much from the results');
end

invE=1:length(E);
invE(E)=invE;

v=(invE(1:n)>invE(n+1:2*n));
