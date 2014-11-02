function [can,invcond]=canBasisFromCanBasis(can,targetp)
% update a canbasis to a new specified "target" p
%
% [can,invcond]=canBasisFromCanBasis(can,targetp)
%
% TODO: update doc equivalent to subspace2CanBasis(canBasis2Subspace(X,p),targetp);
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

[m, n] = size(can.X);

%get the inverse permutation of targetp
invTargetP = 1:n+m;
invTargetP(targetp) = 1:n+m;

%perform a canBasis update that takes the identity submatrix to the correct
%rows (although not necessarily in the correct order)
in = find(invTargetP(can.p(n+1:end))<n+1);
out = find(invTargetP(can.p(1:n))>n);
[can.X can.p invcond]=updateCanBasis(can.X,can.p,in,out);

%reorder the rows and columns of X to get the order right
reordering=invTargetP(can.p);
can.X(reordering(n+1:end)-n,reordering(1:n))=can.X;
