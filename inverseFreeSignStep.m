function [symNew,w,swaps1,swaps2,res,res2]=inverseFreeSignStep(sym,varargin)
% analogous of doublingStep which performs a step of Newton for the matrix sign
%
% [symnew,w,swaps1,swaps2,res,res2]=nmsStep(sym,options)
%
% given X,v symbasis for a Hamiltonian matrix pencil H-sI
% computes H_new=1/2* (H+inv(H))
%
% res and res2 are residual measures. Currently res2 is nan, but it is
% there to get the same signature of doublingStep
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

o=Options(varargin{:});

n=length(sym.X);
[A,E]=hamiltonianPencilFromSymBasis(sym);
Z=[A;E];
[can,invcond1,swaps1]=canBasisFromSubspace(Z,o);
[Etilde,Atilde]=leftDual(can);

%assertVectorsAlmostEqual(Etilde*E,Atilde*A); %I need to be careful with A and E...

%scaling=abs((det(newA)/det(newE)))^(-1/n); %determinantal scaling, see [Higham, function of matrices]
%scaling=norm(newA,'fro')/norm(newE,'fro');

scaling=o.get('scaling',1);
%works with scaling*(E\A) instead of E\A
% the iteration should now read newA=1/2*(scaling*Etilde*A+1/scaling*Atilde*E);
%but this may get us entries which are very large, so we "rescale the
%scaling" and find instead [c,s] such that t*scaling=c, t/scaling=s, and c^2+s^2=1
G=givens(scaling,1/scaling);
s=-G(2,1);c=G(2,2);
assertElementsAlmostEqual(c/s,scaling/inv(scaling));

newA=1/sqrt(2)*(c*Etilde*A+s*Atilde*E);
%newA=1/2*(Etilde*A+Atilde*E);
newE=Etilde*E;

[symNew, invcond2,swaps2]=symBasisFromHamiltonianPencil(newA,newE,o);

%computes residual measures
symOld=symBasisFromSymBasis(sym,symNew.v);
res=norm(symOld.X-symNew.X,'fro');
res2=nan;
w = can.p;
end
