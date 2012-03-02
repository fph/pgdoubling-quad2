function [X,Y,U,V]=solveCARE(A,G,Q,varargin)
% solves a continuous-time algebraic Riccati equation using permuted doubling
%
% [X,Y,U,V]=solveCARE(A,G,Q,...)
%
% options:
% gamma: constant to use for the Cayley transform
% defectCorrectionSteps: makes extra defect correction steps at the end
% (*on X,U only, not Y,V*)
% (default: 0)
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
