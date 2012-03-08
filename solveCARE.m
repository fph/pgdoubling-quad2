function [X,Y,U,V]=solveCARE(A,G,Q,varargin)
% solves a continuous-time algebraic Riccati equation using permuted doubling
%
% [X,Y,U,V]=solveCARE(A,G,Q,...)
%
% options:
% defectCorrectionSteps: makes extra defect correction steps at the end
% (*on X,U only, not Y,V*)
% (default: 0)
% maxSteps, tolerance, verbose: as in doubling.m

o=matgic.Options(varargin{:});

H=hamiltonian(A,G,Q);
%optimizes basis for the Hamiltonian pencil
[S,v]=hamiltonianPencil2SymBasis(H,eye(size(H)));
[S,v]=optimizeSymBasis(S,v);
[AA,EE]=symBasis2HamiltonianPencil(S,v);

gamma=1.1*length(S)*max(max(abs(S))); %this should ensure that gamma does not collide with some eigenvalues 
[S,v]=symplecticPencil2SymBasis(AA+gamma*EE,AA-gamma*EE);
[S,v]=optimizeSymBasis(S,v);

[S,v]=doubling(S,v,'sda',o);

n=length(S)/2;
first=1:n;second=n+1:2*n;

Xpi=-S(second,second);vx=v(second);
Ypi=-S(first,first);vy=v(first);

defectCorrectionSteps=o.get('defectCorrectionSteps',0);

if defectCorrectionSteps>0
    [Xpi,vx]=defectCorrectionCARE(A,G,Q,Xpi,vx,defectCorrectionSteps);
end

U=rowSwap([eye(n);Xpi;],vx,'N');
[X invcond1]=rightLinSolve(U(second,:),U(first,:));

V=rowSwap([Ypi;eye(n)],vy,'T');
[Y invcond2]=rightLinSolve(V(first,:),V(second,:));
