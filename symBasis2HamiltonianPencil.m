function [A,E]=symBasis2HamiltonianPencil(X,v)
% converts a symBasis to a Hamiltonian pencil
%
% [A,E]=symBasis2HamiltonianPencil(X,v)
%
% A,E should be hamiltonian, in that A*J*E
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

U=symBasis2SymplecticSubspace(X,v);
n=length(X);
first=1:n;
second=n+1:2*n;
E=U(first,:)';
A=(jay(n)*U(second,:))';
