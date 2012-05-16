function [X,v,invcond]=symplecticSubspace2SpecifiedSymBasis(U,v)
% symplectic basis representation of a symplectic subspace
%
% [X,v,invcond]=symplecticSubspace2SpecifiedSymBasis(U,v)
%
% given a (matrix whose columns span a) symplectic subspace, returns a symplectic basis representation
% i.e., X and v such that Pi_v * U = [I;X]
% here Pi_v = the matrix in rowSwap(?,v,'N'), and X is symmetric
%
% here, v is specified a priori, no optimization attempt is done
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

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

if invcond<sqrt(eps(class(U))) && nargout<3 %output error only if the caller does not request invcond explicitly. Otherwise, handling it is their business.
    warning('cbrpack:illConditionedMatrix', 'symplecticSubspace2SymBasis: the matrix I am inverting has conditioning >1/sqrt(eps). The selected symplectic swap guess may be wrong');
end
