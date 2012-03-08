function [Xnew,vnew,w,swaps1,swaps2,res,res2]=inverseFreeSignStep(X,v,wguess,vguess,threshold1,threshold2d,threshold2o)
% analogous of doublingStep which performs a step of Newton for the matrix sign
%
% [Xnew,vnew,w,swaps1,swaps2,res,res2]=nmsStep(X,v,wguess,vguess,threshold1,threshold2d,threshold2o)
%
% given X,v symbasis for a Hamiltonian matrix pencil H-sI
% computes H_new=1/2* (H+inv(H))
%
% res and res2 are residual measures. Currently res2 is nan, but it is
% there to get the same signature of doublingStep

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
[leftX,w,swaps1,optcond]=optimizeCanBasis(leftX,w,threshold1);
if 1/optcond>n*threshold1
    %'recompute' --- this almost never happens
    [leftX,w,invcond1]=subspace2CanBasis(Z,w);
end
[Etilde,Atilde]=leftDual(leftX,w);

%assertVectorsAlmostEqual(Etilde*E,Atilde*A); %I need to be careful with A and E...

newA=1/2*(Etilde*A+Atilde*E);
newE=Etilde*E;

%scaling=abs((det(newA)/det(newE)))^(-1/n); %determinantal scaling, see [Higham, function of matrices]
%scaling=norm(newA,'fro')/norm(newE,'fro');
% this is an odd place to scale. Usually you want to scale directly in the
% NMS iteration, H<- 1/2(s*H+inv(s*H))
% But when you go inverse-free, this does not work as it makes the matrices
% unbounded (s may be very large/small). Instead, we scale here, which is
% equivalent but plays better with the inverse-free setting.
%
% TODO: disabled scaling because it gave problems. Need to investigate
[Xnew vnew invcond2]=hamiltonianPencil2SymBasis(newA,newE,vguess);
[Xnew vnew swaps2,optcond]=optimizeSymBasis(Xnew,vnew,threshold2d,threshold2o);
if 1/optcond>n*threshold2o
    %'recompute' --- this almost never happens
    [Xnew vnew invcond2]=hamiltonianPencil2SymBasis(newA,newE,vnew);
end
%computes residual measures
Xold=symBasis2symBasis(X,v,vnew);
res=norm(Xold-Xnew,'fro');
res2=nan;
end
