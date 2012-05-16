function [X,p,invcond]=subspace2SpecifiedCanBasis(U,p)
% compute a fixed canonical basis representation of a subspace
%
% [X,p,invcond]=subspace2CanBasis(U,p);
%
% input: U subspace, p permutation to use
%
% output: 
% X,p: such that U(p,:) and [I;X] span the same subspace
% p is copied to the output without changes (it is there to ensure
% consistent output with the other procedures)
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

[m n]=size(U);
assertEqual(length(p),m);

%X=U(p(n+1:m),:)/U(p(1:n),:);
[X,invcond]=rightLinSolve(U(p(n+1:m),:),U(p(1:n),:));

if invcond<sqrt(eps(class(U))) && nargout<3
    warning('cbrpack:illConditionedMatrix', 'subspace2CanBasis: the matrix I am inverting has conditioning >1/sqrt(eps). This may be due to an ill-conditioned subspace or to a bad initial guess --- consider using the initial value heuristic instead');
end
