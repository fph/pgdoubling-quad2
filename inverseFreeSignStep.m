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
[A,E]=symBasis2HamiltonianPencil(X,v);
Z=[A;E];
[leftX,w,invcond1,swaps1]=subspace2CanBasis(Z,'threshold',threshold1,'initialPermutation',wguess);
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
[Xnew vnew invcond2,swaps2]=hamiltonianPencil2SymBasis(newA,newE,'diagonalThreshold',threshold2d,'offDiagonalThreshold',threshold2o,'initialSwap',vguess);

%computes residual measures
Xold=symBasis2symBasis(X,v,vnew);
res=norm(Xold-Xnew,'fro');
res2=nan;
end
