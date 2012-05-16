function H=hamiltonian(A,G,Q)
% builds the Hamiltonian of an ARE
%
% H=hamiltonian(A,G,Q)
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

H=[A -G;-Q -A'];
