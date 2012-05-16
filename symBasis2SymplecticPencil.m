function [L,U]=symBasis2SymplecticPencil(X,v)
% unpacks the symBasis representation of a symplectic pencil
%
% function [L,U]=symBasis2SymplecticPencil(X,v)
%
% given X,v representing a symplectic pencil,
% returns L,U such that pencil=L-s*U
% L is permuted lower block triangular, U upper
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

n=length(X)/2;
%if we ever adapt this to a nonsymmetric case, add a X=X' here
first=1:n;second=n+1:2*n;

U=[[eye(n);zeros(n)] X(:,first)];

L=[X(:,second) [zeros(n);eye(n)]];

%undo the column swaps
U=rowSwap(U',v(first),'T');
U=U';

L=rowSwap(L',v(second),'N');
L=L';
