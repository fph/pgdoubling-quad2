function [X,Y,U,V]=solveEcareSign(A,B,Q,R,S,varargin)
% solves a CARE given in 3x3 even block form using inverse-free sign
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


[S,v]=inverseFreeSign(S,v,o);

[A,E]=symBasis2HamiltonianPencil(S,v);

%to be replaced by sth else structure-preserving...
n=length(S)/2;
first=1:n;second=n+1:2*n;
[u s v]=svd(A+E);U=v(:,second);
[u s v]=svd(A-E);V=v(:,second);


[X invcond1]=rightLinSolve(U(first,:),U(second,:));
[Y invcond2]=rightLinSolve(V(second,:),V(first,:));
