function [A,E]=hamiltonianPencilFromSymBasis(sym)
% converts a symBasis to a Hamiltonian pencil
%
% [A,E]=hamiltonianPencilFromSymBasis(sym)
%
% The returned pencil is Hamiltonian.
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling


if not(strcmp(sym.origin,'hamiltonianPencil'))
    warning('PGDoubling:wrongOrigin','You are converting to a Hamiltonian pencil a symBasis that did not come from one.');
end
sym.origin='symplecticSubspace';
U=symplecticSubspaceFromSymBasis(sym);
n=length(sym.X);
first=1:n;
second=n+1:2*n;
E=U(first,:)';
A=(jay(n)*U(second,:))';
