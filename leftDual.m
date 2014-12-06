function [Etilde,Atilde]=leftDual(can)
% computes left dual of a subspace (given its canBasis)
%
% [Etilde,Atilde]=leftDual(can)
%
% given a canBasis of U (of dim 2nxn), computes Etilde,Atilde such that 
% [Etilde,Atilde]*J*U=0 and rk [Etilde,Atilde]=n
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

n=length(can.X);
first=1:n;second=n+1:2*n;
can.X=-can.X';
can.p=can.p([second,first]);
Tildes=subspaceFromCanBasis(can);Tildes=Tildes';
Etilde=Tildes(:,second);Atilde=-Tildes(:,first);
