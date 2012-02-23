function [X,Y,U,V]=solveECARE(A,B,Q,R,S,varargin)
% solves a CARE given in 3x3 even block form, and computes invariant subspaces
%
% [X,Y,U,V]=solveECARE(A,B,Q,R,S,varargin)
%
% returns the first two blocks of the stable and unstable invariant subspaces 
% of evenPencil(A,B,C,Q,R,S), U and V.
%
% X and Y are obtained by writing them in form [X;I] and [I;Y]
% respectively. In other words, X is the stabilizing solution of the Riccati equation
% with coefficients A-BR^(-1)C', -BR^(-1)B', Q-CR^(-1)C', and Y the one of
% its dual equation.
%
% options as in solveCARE

o=matgic.Options(varargin{:});

[n m]=size(B);

if not(exist('S','var'))
    S=[];
end

if not(exist('v','var'))
    v=[];
end

[AA,EE]=evenPencil(A,B,Q,R,S);

[S,v]=evenPencil2SymBasis(AA,EE,n,m,v);
[S,v]=optimizeSymBasis(S,v);

[Ah,Eh]=symBasis2HamiltonianPencil(S,v);

gamma=o.get('gamma',norm(Ah,'fro')/norm(Eh,'fro'));
if not(gamma>0)
    error 'gamma must be positive'
end

[S,v]=symplecticPencil2SymBasis(Ah+gamma*Eh,Ah-gamma*Eh);
[S,v]=optimizeSymBasis(S,v);
%TODO: need to swap something to adjust [I;X] vs. [X;I]
[X,Y,U,V]=doubling(S,v,o);

