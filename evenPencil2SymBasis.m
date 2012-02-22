function [X,v]=evenPencil2SymBasis(AA,EE,n,m,v)
%converts an even pencil to a symBasis
%
% [X,v]=evenPencil2SymBasis(AA,EE,n,m,v)
%
% since we wish to allow more general pencils AA-sEE, we give n and m
% explicitly
%
% Warning: kinda misnamed, does not work with every even pencil --- I need
% to write down the theory

firstSecond=1:2*n;third=2*n+1:2*n+m;
U=[EE(:,firstSecond)'; (AA(:,firstSecond)*jay(2*n))'];

U=rowSwap(U,v,'N');

Xext=U(2*n+1:4*n,:)/[U(1:2*n,:); AA(third,:)];

X=Xext(:,firstSecond);
