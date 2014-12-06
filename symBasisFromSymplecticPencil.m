function [sym,invcond,swaps]=symBasisFromSymplecticPencil(L,U,varargin)
% constructs a canonical basis representation of a symplectic pencil
%
% [sym,swaps]=symBasisFromSymplecticPencil(L,U,varargin)
%
% L,U define pencil L-sU which is symplectic, i.e., L*jay(length(L))*L'=U*jay(length(U))*U',
% optional arguments and outputs are as in symplecticSubspace2SymBasis.m
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

o=Options(varargin{:});

n=length(L)/2;
first=1:n;second=n+1:2*n;

[sym,invcond,swaps]=symBasisFromSymplecticSubspace([U(:,first) L(:,second) U(:,second) L(:,first)]',o);
sym.origin='symplecticPencil';
