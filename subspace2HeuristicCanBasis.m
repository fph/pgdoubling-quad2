function [X,p,invcond]=subspace2HeuristicCanBasis(U)
% compute a fixed canonical basis representation of a subspace
%
% [X,p,invcond]=subspace2HeuristicCanBasis(U);
%
% input: U subspace, p permutation to use
%
% output: 
% X,p: such that U(p,:) and [I;X] span the same subspace
%
% p is determined through QRP, and the same QRP is used to solve the system

[m n]=size(U);

[Q,R,p]=qr(U',0);

R=R';

[X,invcond]=rightLinSolve(R(n+1:m,:),R(1:n,:));

if invcond<sqrt(eps(class(U))) && nargout<3
    warning('cbrpack:illConditionedMatrix', 'subspace2CanBasis: the matrix I am inverting has conditioning >1/sqrt(eps). This may be due to an ill-conditioned subspace or to a bad initial guess --- consider using the initial value heuristic instead');
end
