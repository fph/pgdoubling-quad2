function [Xnew,vnew,w,swaps1,swaps2,res]=inverseFreeSignStep(X,v,wguess,vguess,threshold1,threshold2d,threshold2o)
% analogous of doublingStep which performs a step of Newton for the matrix sign
%
% [Xnew,vnew,w,swaps1,swaps2,res]=nmsStep(X,v,wguess,vguess,threshold1,threshold2d,threshold2o)
%
% given X,v symbasis for a Hamiltonian matrix pencil H-sI
% computes H_new=1/2* (H+inv(H))
%

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
first=1:n;second=n+1:2*n;
[A,E]=symBasis2HamiltonianPencil(X,v);
Z=[A;E];
[leftX,w,invcond1]=subspace2CanBasis(Z,wguess);
[leftX,w,swaps1]=optimizeCanBasis(leftX,w,threshold1);
[leftX,w,invcond1]=subspace2CanBasis(Z,w);
[Etilde,Atilde]=leftDual(leftX,w);
assertVectorsAlmostEqual(Etilde*E,Atilde*A); %I need to be careful with A and E...
newA=1/2*(Etilde*A+Atilde*E);
newE=Etilde*E;
%TODO: scaling
[Xnew vnew invcond2]=hamiltonianPencil2SymBasis(newA,newE,vguess);
[Xnew vnew swaps2]=optimizeSymBasis(Xnew,vnew,threshold2d,threshold2o);
[Xnew vnew invcond2]=hamiltonianPencil2SymBasis(newA,newE,vnew);

%computes residual measures
Xold=symBasis2symBasis(X,v,vnew);
res=norm(Xold-Xnew,'fro');
end
