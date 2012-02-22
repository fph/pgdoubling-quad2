function [A,E]=symBasis2HamiltonianPencil(X,v)
% converts a symBasis to a Hamiltonian pencil
%
% [A,E]=symBasis2HamiltonianPencil(X,v)
%
% A,E should be hamiltonian, in that A*J*E

U=symBasis2SymplecticSubspace(X,v);
n=length(X);
first=1:n;
second=n+1:2*n;
E=U(first,:)';
A=(jay(n)*U(second,:))';
