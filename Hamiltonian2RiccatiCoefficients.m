function [A,G,Q]=Hamiltonian2RiccatiCoefficients(H)
% converts a Hamiltonian H to Riccati coefficients
%
% [A,G,Q]=Hamiltonian2RiccatiCoefficients(H)

n=length(H)/2;
first=1:n;second=n+1:2*n;
A=H(first,first);
G=-H(first,second);
Q=-H(second,first);
