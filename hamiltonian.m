function H=hamiltonian(A,G,Q)
% builds the Hamiltonian of an ARE
%
% H=hamiltonian(A,G,Q)

H=[A -G;-Q -A'];
