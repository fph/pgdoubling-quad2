function [can,invcond]=specifiedCanBasisFromSubspace(U,p)
% compute a fixed canonical basis representation of a subspace
%
% [can,invcond]=specifiedCanBasisFromSubspace(U,p);
%
% input: U subspace, p permutation to use
%
% output: 
% can: such can.p == p and U(can.p,:) and [I;can.X] span the same subspace
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

[m n]=size(U);
if m ~= length(p)
    error('PGDoubling:wrongDimension', 'wrong length of the permutation vector');
end

%X=U(p(n+1:m),:)/U(p(1:n),:);
[can.X,invcond] = rightLinSolve(U(p(n+1:m),:),U(p(1:n),:));
can.p = p;
can.origin = 'subspace';

if invcond<sqrt(eps(class(U))) && nargout<3
    warning('PGDoubling:illConditionedMatrix', 'CanBasisFromSubspace: the matrix I am inverting has conditioning >1/sqrt(eps). This may be due to an ill-conditioned subspace or to a bad initial guess --- consider using the initial value heuristic instead');
end
