function [X,v,swaps]=symplecticPencil2SymBasis(L,U,varargin)
% constructs a canonical basis representation of a symplectic pencil
%
% [X,v,swaps]=symplecticPencil2SymBasis(L,U,...)
%
% L,U define a symplectic pencil, i.e. L*jay(length(L))*L'=U*jay(length(U))*U',
% optional arguments and outputs are as in symplecticSubspace2SymBasis.m

n=length(L)/2;
first=1:n;second=n+1:2*n;

[X,v,swaps]=symplecticSubspace2SymBasis([U(:,first) L(:,second) U(:,second) L(:,first)]',varargin{:});
