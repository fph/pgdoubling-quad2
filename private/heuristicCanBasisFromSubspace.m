function [can,invcond]=heuristicCanBasisFromSubspace(U)
% compute a fixed canonical basis representation of a subspace
%
% [can,invcond]=heuristicCanBasisFromSubspace(U);
%
% input: U subspace
%
% output: 
% can: such that U(can.p,:) and [I;can.X] span the same subspace
% invcond: inverse condition number
%
% p is determined through QRP, and the same QRP is used to solve the system
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

[m, n] = size(U);

[Q, R, can.p] = qr(U', 0);

R = R';

[can.X, invcond] = rightLinSolve(R(n+1:m,:),R(1:n,:));
can.origin = 'subspace';

if invcond < sqrt(eps(class(U))) && nargout < 3
    warning('cbrpack:illConditionedMatrix2', 'subspace2CanBasis: the matrix I am inverting has conditioning >1/sqrt(eps). This is probably due to an ill-conditioned subspace');
end
