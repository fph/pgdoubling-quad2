function [Xnew,vnew,swaps1,swaps2,w,nn]=doublingStep(X,v,options1,options2)
% perform a doubling step using permuted bases
%
% [Xnew,vnew,swaps1,swaps2,w,nn]=doublingStep(X,v,options1,options2)
%
%
% (X,v) is a symBasis for a symplectic pencil L-sU
% options1 are the options to use for pencil2canBasis(L,U) 
% options2 are the options to use for symplecticpencil2canBasis(Lsquare,Usquare)
% seriously consider specifying initial guesses when available
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

if not(exist('options1','var')) || isempty(options1)
    options1=matgic.Options();
end

if not(exist('options2','var')) || isempty(options2)
    options2=matgic.Options();
end

n=length(X);
first=1:n;second=n+1:2*n;
[L,U]=symBasis2SymplecticPencil(X,v);
Z=[L;U];
[leftX,w,swaps1]=subspace2CanBasis(Z,options1);
invw=w([second,first]);
Tildes=canBasis2Subspace(-leftX',invw);Tildes=Tildes';
Ltilde=Tildes(:,second);Utilde=-Tildes(:,first);
assertVectorsAlmostEqual(Ltilde*U,Utilde*L);
newL=Ltilde*L;
newU=Utilde*U;
%TODO: could do scaling as in Newton 
[Xnew vnew swaps2]=symplecticPencil2SymBasis(newL,newU,options2);
if all(vnew(:)==v(:))
    n=n/2;first=1:n;second=n+1:2*n;
    nn=norm(S(first,first)-S_new(first,first),'fro')+norm(S(second,second)-S_new(second,second),'fro');
else
    nn=nan;
end
