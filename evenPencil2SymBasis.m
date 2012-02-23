function [X,v,invcond]=evenPencil2SymBasis(AA,EE,n,m,v)
%converts an even pencil to a symBasis
%
% [X,v,invcond]=evenPencil2SymBasis(AA,EE,n,m,v)
%
% since we wish to allow more general pencils AA-sEE, we give n and m
% explicitly
%
% Warning: kinda misnamed, does not work with every even pencil --- I need
% to write down the theory

firstSecond=1:2*n;third=2*n+1:2*n+m;
U=[EE(:,firstSecond)'; (AA(:,firstSecond)*jay(2*n))'];

if not(exist('v','var')) || isempty(v)
    v=extendedSymBasisHeuristic([U; AA(:,third)']);
end

U=rowSwap(U,v,'N');


[Xext invcond]=rightLinSolve(U(2*n+1:4*n,:),[U(1:2*n,:); AA(third,:)]);

X=Xext(:,firstSecond);

%final check
if norm(X-X','fro')/norm(X,'fro') > sqrt(eps)
    warning('cbrpack:notSymplectic','the resulting matrix is numerically very far from Hermitian --- was your starting subspace symplectic?');
end

%symmetrize
X=(X+X')/2;

if(invcond<sqrt(eps(class(U))))
    warning('cbrpack:illConditionedMatrix', 'symplecticSubspace2SymBasis: the matrix I am inverting has conditioning >1/sqrt(eps). This may be due to an ill-conditioned subspace or to a bad initial guess --- consider using the initial value heuristic instead');
end
