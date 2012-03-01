function [X,Y,U,V]=solveECARE(A,B,Q,R,S,varargin)
% solves a CARE given in 3x3 even block form, and computes invariant subspaces
%
% [X,Y,U,V]=solveECARE(A,B,Q,R,S,varargin)
%
% returns the first two blocks of the stable and unstable invariant subspaces 
% of evenPencil(A,B,C,Q,R,S), U and V.
%
% Accepts options as in solveCARE
%
% X and Y are obtained by writing them in form [X;I] and [I;Y]
% respectively. In other words, X is the stabilizing solution of the Riccati equation
% with coefficients A-BR^(-1)C', -BR^(-1)B', Q-CR^(-1)C', and Y the one of
% its dual equation.
%
% Notice though that the order of the blocks in the ECARE and the CARE
% approach are different, so U and V obtained by this function and
% solveCARE differ by a factor abs(jay()), and we get [X;I] instead of
% [I;X]. This is a fault of the standard notations, not mine.

o=matgic.Options(varargin{:});

[n m]=size(B);

if not(exist('S','var'))
    S=[];
end

v=o.get('initialv',[]);

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

[S,v]=doubling(S,v,o);

n=length(S)/2;
first=1:n;second=n+1:2*n;

U=rowSwap([eye(n);-S(second,second);],v(second),'N');
[X invcond1]=rightLinSolve(U(first,:),U(second,:));

V=rowSwap([-S(first,first);eye(n)],v(first),'T');
[Y invcond2]=rightLinSolve(V(second,:),V(first,:));