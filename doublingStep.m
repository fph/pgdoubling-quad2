function [Xnew,vnew,w,swaps1,swaps2,nGH,nEF]=doublingStep(X,v,wguess,vguess,threshold1,threshold2d,threshold2o)
% perform a doubling step using permuted bases
%
% [Xnew,vnew,swaps1,swaps2,w,nn]=doublingStep(X,v,wguess,vguess,threshold1,threshold2d,threshold2o)
%
% (X,v) is a symBasis for a symplectic pencil L-sU
% 
% (Xnew,vnew,swaps2) (in output) is a symBasis for the symplectic pencil obtained by "squaring" L-sU
% squaring here is extended to pencils in the sense of [Benner, Byers]:
% if U is invertible, then it is a pencil that is right-equivalent to U^{-1}LU^{-1}L-sI,
% and this is extended by continuity to singular U.
% 
% w,swaps1 are the permutation and swap resulting from pencil2CanBasis(A,E)
% (w can be recycled in options1 for a successive step)
% 
% note that one swap in swaps2 should cost half as much as one in swaps1, because of symmetry.
%
% if vnew==v, also computes a cheap measure of the change in the subspace that can be used as a stopping criterion

if not(exist('wguess','var'))
    wguess=[];
end

if not(exist('vguess','var'))
    vguess=[];
end

if not(exist('threshold1','var'))
    threshold1=[];
end
if not(exist('threshold2d','var'))
    threshold2d=[];
end
if not(exist('threshold2o','var'))
    threshold2o=[];
end
 
n=length(X);
[L,U]=symBasis2SymplecticPencil(X,v);
Z=[L;U];
[leftX,w,invcond1]=subspace2CanBasis(Z,wguess);
[leftX,w,swaps1,optcond]=optimizeCanBasis(leftX,w,threshold1);
if 1/optcond>n*threshold1
    %'recompute' --- this almost never happens
    [leftX,w,invcond1]=subspace2CanBasis(Z,w);
end

[Ltilde,Utilde]=leftDual(leftX,w);
%assertVectorsAlmostEqual(Ltilde*U,Utilde*L);
newL=Ltilde*L;
newU=Utilde*U;
%TODO: could do scaling as in Newton for the matrix sign

[Xnew vnew invcond2]=symplecticPencil2SymBasis(newL,newU,vguess);
[Xnew vnew swaps2 optcond]=optimizeSymBasis(Xnew,vnew,threshold2d,threshold2o);
if 1/optcond>n*threshold2o
    %'recompute' --- this almost never happens
    [Xnew vnew invcond2]=symplecticPencil2SymBasis(newL,newU,vnew);
end

%computes residual measures
Xold=symBasis2symBasis(X,v,vnew);
n=n/2;first=1:n;second=n+1:2*n;
nGH=norm(Xold(first,first)-Xnew(first,first),'fro')+norm(Xold(second,second)-Xnew(second,second),'fro');
nEF=norm(Xnew(first,second),'fro')+norm(Xnew(second,first),'fro');
end
