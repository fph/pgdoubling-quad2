function [Etilde,Atilde]=leftDual(X,p)
% computes left dual of a subspace (given its canBasis)
%
% [Etilde,Atilde]=leftDual(X,p)
%
% given a canBasis of U (of dim 2nxn), computes Etilde,Atilde such that 
% [Etilde,Atilde]*J*U=0 and rk [Etilde,Atilde]=n

n=length(X);
first=1:n;second=n+1:2*n;
invp=p([second,first]);
Tildes=canBasis2Subspace(-X',invp);Tildes=Tildes';
Etilde=Tildes(:,second);Atilde=-Tildes(:,first);
