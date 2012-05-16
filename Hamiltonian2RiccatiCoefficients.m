function [A,G,Q]=Hamiltonian2RiccatiCoefficients(H)
% converts a Hamiltonian H to Riccati coefficients
%
% [A,G,Q]=Hamiltonian2RiccatiCoefficients(H)
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

n=length(H)/2;
first=1:n;second=n+1:2*n;
A=H(first,first);
G=-H(first,second);
Q=-H(second,first);
