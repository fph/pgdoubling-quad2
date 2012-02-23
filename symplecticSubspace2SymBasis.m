function [X,v,invcond]=symplecticSubspace2SymBasis(U,v)
% symplectic basis representation of a symplectic subspace
%
% [X,v,invcond]=symplecticSubspace2SymBasis(U,v)
%
% given a (matrix whose columns span a) symplectic subspace, returns a symplectic basis representation
% i.e., X and v such that Pi_v * U = [I;X]
% here Pi_v = the matrix in rowSwap(?,v,'N'), and X is symmetric
%
%
% The basis given may not meet a given threshold --- run optimizeSymBasis
% if you need an optimal basis
%
% If your initial guess is off, you might need to do this twice to get more
% stability:
%
% [X,v]=symplecticSubspace2SymBasis(U,crappyInitialGuess);
% [X,v]=optimizeSymBasis(X,v);
% [X,v]=symplecticsubspace2SymBasis(U,v); %should be ok now!
% 

if not(exist('v','var')) || isempty(v)
    v=symBasisHeuristic(U);
end

[m n]=size(U);
if m~=2*n
    error('cbrpack:oddSize','symplectic matrices should be 2nxn');
end

S=rowSwap(U,v,'N');

[X invcond]=rightLinSolve(S(n+1:end,:),S(1:n,:));

%final check
if norm(X-X','fro')/norm(X,'fro') > sqrt(eps)
    warning('cbrpack:notSymplectic','the resulting matrix is numerically very far from Hermitian --- was your starting subspace symplectic?');
end

%symmetrize
X=(X+X')/2;

if(invcond<sqrt(eps(class(U))))
    warning('cbrpack:illConditionedMatrix', 'symplecticSubspace2SymBasis: the matrix I am inverting has conditioning >1/sqrt(eps). This may be due to an ill-conditioned subspace or to a bad initial guess --- consider using the initial value heuristic instead');
end
