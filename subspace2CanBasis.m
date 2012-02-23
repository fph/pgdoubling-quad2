function [X,p,invcond]=subspace2CanBasis(U,p)
% Canonical basis representation of a subspace
%
% [X,p,invcond]=subspace2CanBasis(U,p);
%
% input: U subspace, p initial permutation guess (may be empty, in this 
% case the procedure will use a O(n^3) heuristic based on QRP)
%
% output: 
% X,p: such that U(p(1:end),:) and [I;X] span the same subspace
%
% The basis given may not meet a given threshold --- run optimizeCanBasis
% if you need an optimal basis
%
% If your initial guess is off, you might need to do this twice to get more
% stability:
%
% [X,p]=subspace2CanBasis(U,crappyInitialGuess);
% [X,p]=optimizeCanBasis(X,p);
% [X,p]=subspace2CanBasis(U,p); %should be ok now!

if not(exist('p','var')) || isempty(p)
    p=canBasisHeuristic(U);
end

[m n]=size(U);
assertEqual(length(p),m);

%X=U(p(n+1:m),:)/U(p(1:n),:);
[X,invcond]=rightLinSolve(U(p(n+1:m),:),U(p(1:n),:));

if(invcond<sqrt(eps(class(U))))
    warning('cbrpack:illConditionedMatrix', 'subspace2CanBasis: the matrix I am inverting has conditioning >1/sqrt(eps). This may be due to an ill-conditioned subspace or to a bad initial guess --- consider using the initial value heuristic instead');
end
