function [X,Y,U,V]=solveCARE(A,G,Q,varargin)
% solves a continuous-time algebraic Riccati equation using permuted doubling
%
% [X,Y,U,V]=solveCARE(A,G,Q,...)
%
% options:
% gamma: constant to use for the Cayley transform
% maxSteps, tolerance, verbose: as in doubling.m

o=matgic.Options(varargin{:});

H=hamiltonian(A,G,Q);
gamma=o.get('gamma',norm(H));
if not(gamma>0)
    error 'gamma must be positive'
end
[S,v]=symplecticPencil2SymBasis(H+gamma*eye(size(H)),H-gamma*eye(size(H)));
[S,v]=optimizeSymBasis(S,v);

[S,v]=doubling(S,v,o);

n=length(S)/2;
first=1:n;second=n+1:2*n;

U=rowSwap([eye(n);-S(second,second);],v(second),'N');
[X invcond1]=rightLinSolve(U(second,:),U(first,:));

V=rowSwap([-S(first,first);eye(n)],v(first),'T');
[Y invcond2]=rightLinSolve(V(first,:),V(second,:));
