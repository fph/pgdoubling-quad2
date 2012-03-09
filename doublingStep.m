function [Xnew,vnew,w,swaps1,swaps2,res,res2]=doublingStep(X,v,varargin)
% perform a doubling step using permuted bases
%
% [Xnew,vnew,w,swaps1,swaps2,res,res2]=doublingStep(X,v,options)
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
% also computes two residual measures

o=matgic.Options(varargin{:});
 
n=length(X);
[L,U]=symBasis2SymplecticPencil(X,v);
Z=[L;U];
[leftX,w,invcond1,swaps1]=subspace2CanBasis(Z,o);

[Ltilde,Utilde]=leftDual(leftX,w);
%assertVectorsAlmostEqual(Ltilde*U,Utilde*L);
newL=Ltilde*L;
newU=Utilde*U;
%TODO: could do scaling as in Newton for the matrix sign

[Xnew vnew invcond2 swaps2]=symplecticPencil2SymBasis(newL,newU,o);

%computes residual measures
Xold=symBasis2symBasis(X,v,vnew);
res2=norm(Xold-Xnew,'fro');
n=n/2;first=1:n;second=n+1:2*n;
res=norm(Xnew(first,second),'fro')+norm(Xnew(second,first),'fro');
end
