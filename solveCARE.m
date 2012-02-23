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

[X,Y,U,V]=doubling(S,v,o);
