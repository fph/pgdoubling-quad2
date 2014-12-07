function [sym2,w,swaps1,swaps2,res,res2]=doublingStep(sym,varargin)
% perform a doubling step using permuted bases
%
% [sym2,w,swaps1,swaps2,res,res2]=doublingStep(sym,options)
%
% (X,v) is a symBasis for a symplectic pencil L-sU
% 
% (sym2,swaps2) (in output) is a symBasis for the symplectic pencil obtained by "squaring" L-sU
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
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

o=Options(varargin{:});
 
n=length(sym.X);
[L,U]=symplecticPencilFromSymBasis(sym);
Z=[L;U];
[can,invcond1,swaps1]=canBasisFromSubspace(Z,o);

[Ltilde,Utilde]=leftDual(can);
%assertVectorsAlmostEqual(Ltilde*U,Utilde*L);
newL=Ltilde*L;
newU=Utilde*U;
%it is not easy to do scaling as in Newton for the matrix sign, in this
%setting

[sym2 invcond2 swaps2]=symBasisFromSymplecticPencil(newL,newU,o);

%computes residual measures
symOld=symBasisFromSymBasis(sym,sym2.v);
res2=norm(symOld.X-sym2.X,'fro');
n=n/2;first=1:n;second=n+1:2*n;
res=norm(sym2.X(first,second),'fro')+norm(sym2.X(second,first),'fro');
w = can.p;
end
