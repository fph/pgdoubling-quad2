function [X,v,invcond]=symplecticSubspace2HeuristicSymBasis(U)
% symplectic basis representation of a symplectic subspace
%
% [X,v,invcond]=symplecticSubspace2SpecifiedSymBasis(U)
%
% given a (matrix whose columns span a) symplectic subspace, returns a symplectic basis representation
% i.e., X and v such that Pi_v * U = [I;X]
% here Pi_v = the matrix in rowSwap(?,v,'N'), and X is symmetric
%
% v is obtained by the PiRQ heuristic

[m n]=size(U);
if m~=2*n
    error('cbrpack:oddSize','symplectic matrices should be 2nxn');
end

[v,R,Q,invcond]=pirq(U);
[X invcond]=rightLinSolve(R(n+1:end,:),R(1:n,:));

%final check
if norm(X-X','fro')/norm(X,'fro') > sqrt(eps)
    warning('cbrpack:notSymplectic','the resulting matrix is numerically very far from Hermitian --- was your starting subspace symplectic?');
end

%symmetrize
X=(X+X')/2;

if invcond<sqrt(eps(class(U))) && nargout<3 %output error only if the caller does not request invcond explicitly. Otherwise, handling it is their business.
    warning('cbrpack:illConditionedMatrix', 'symplecticSubspace2SymBasis: the matrix I am inverting has conditioning >1/sqrt(eps). The selected symplectic swap guess may be wrong');
end
