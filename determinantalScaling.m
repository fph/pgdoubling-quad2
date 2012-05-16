function mu=determinantalScaling(S,v);
% returns the determinantal scaling factor for the matrix sign iteration
%
% mu=determinantalScaling(S,v);
%
% returns mu=abs(det(E)/det(A))^(1/n), where [A,E] is the Hamiltonian
% pencil with symBasis S,v
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

[AA,EE]=symBasis2HamiltonianPencil(S,v);
[L,U]=lu(AA);
aa=exp(mean(log(abs(diag(U)))));

[L,U]=lu(EE);
ee=exp(mean(log(abs(diag(U)))));

mu=ee/aa;

