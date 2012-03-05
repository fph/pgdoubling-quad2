function [X,Y,U,V]=solveCareSign(A,G,Q,varargin)
% solves a CARE with inverse-free sign
%
% [X,Y,U,V]=solveCareSign(A,G,Q,varargin);

o=matgic.Options(varargin{:});

H=hamiltonian(A,G,Q);
[S,v]=hamiltonianPencil2SymBasis(H,eye(size(H)));

[S,v]=inverseFreeSign(S,v,o);

[A,E]=symBasis2HamiltonianPencil(S,v);

%to be replaced by sth else structure-preserving...
n=size(S)/2;
first=1:n;second=n+1:2*n;
[u s v]=svd(A+E);U=v(:,second);
[u s v]=svd(A-E);V=v(:,second);


[X invcond1]=rightLinSolve(U(second,:),U(first,:));
[Y invcond2]=rightLinSolve(V(first,:),V(second,:));
