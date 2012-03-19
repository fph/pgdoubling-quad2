function mu=determinantalScaling(S,v);
% returns the determinantal scaling factor for the matrix sign iteration
%
% mu=determinantalScaling(S,v);
%
% returns mu=abs(det(E)/det(A))^(1/n), where [A,E] is the Hamiltonian
% pencil with symBasis S,v

[AA,EE]=symBasis2HamiltonianPencil(S,v);
[L,U]=lu(AA);
aa=exp(mean(log(abs(diag(U)))));

[L,U]=lu(EE);
ee=exp(mean(log(abs(diag(U)))));

mu=ee/aa;

