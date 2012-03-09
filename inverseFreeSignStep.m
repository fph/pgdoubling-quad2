function [Xnew,vnew,w,swaps1,swaps2,res,res2]=inverseFreeSignStep(X,v,varargin)
% analogous of doublingStep which performs a step of Newton for the matrix sign
%
% [Xnew,vnew,w,swaps1,swaps2,res,res2]=nmsStep(X,v,options)
%
% given X,v symbasis for a Hamiltonian matrix pencil H-sI
% computes H_new=1/2* (H+inv(H))
%
% res and res2 are residual measures. Currently res2 is nan, but it is
% there to get the same signature of doublingStep

o=matgic.Options(varargin{:});

n=length(X);
[A,E]=symBasis2HamiltonianPencil(X,v);
Z=[A;E];
[leftX,w,invcond1,swaps1]=subspace2CanBasis(Z,o);
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
[Xnew vnew invcond2,swaps2]=hamiltonianPencil2SymBasis(newA,newE,o);

%computes residual measures
Xold=symBasis2symBasis(X,v,vnew);
res=norm(Xold-Xnew,'fro');
res2=nan;
end
