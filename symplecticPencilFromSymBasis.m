function [L,U]=symplecticPencilFromSymBasis(sym)
% unpacks the symBasis representation of a symplectic pencil
%
% function [L,U]=symplecticPencilFromSymBasis(sym)
%
% given sym representing a symplectic pencil,
% returns L,U such that pencil=L-s*U
% L is permuted lower block triangular, U upper
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

if not(strcmp(sym.origin,'symplecticPencil'))
    warning('PGDoubling:wrongOrigin','You are converting to a symplectic pencil a symbasis that did not originate from a symplectic pencil');
end

n=length(sym.X)/2;
%if we ever adapt this to a nonsymmetric case, add a X=X' here
first=1:n;second=n+1:2*n;

U=[[eye(n);zeros(n)] sym.X(:,first)];

L=[sym.X(:,second) [zeros(n);eye(n)]];

%undo the column swaps
U=rowSwap(U',sym.v(first),'T');
U=U';

L=rowSwap(L',sym.v(second),'N');
L=L';
