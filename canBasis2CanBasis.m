function [X,targetp,invcond]=canBasis2CanBasis(X,p,targetp)
% update a canbasis to a new specified "target" p
%
% [X,targetp,invcond]=canBasis2CanBasis(X,p,targetp)
%
% equivalent to subspace2CanBasis(canBasis2Subspace(X,p),targetp);
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

[m n]=size(X);

%get the inverse permutation of targetp
invTargetP=1:n+m;
invTargetP(targetp)=1:n+m;

%perform a canBasis update that takes the identity submatrix to the correct
%rows (although not necessarily in the correct order)
in=find(invTargetP(p(n+1:end))<n+1);
out=find(invTargetP(p(1:n))>n);
[X pnew invcond]=updateCanBasis(X,p,in,out);

%reorder the rows and columns of X to get the order right
reordering=invTargetP(pnew);
X(reordering(n+1:end)-n,reordering(1:n))=X;
