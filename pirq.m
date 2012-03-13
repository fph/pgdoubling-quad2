function [v,U,Q,invcond]=pirq(U)
%computes a Pi*R*Q factorization of a symplectic subspace U
%
% [v,U,Q,invcond]=pirq(U)
%
% U has size 2n x n
%
% returns Pi_v, R, Q such that R*Q=RowSwap(U,v,'N')
% R is permuted lower triangular, Q is orthogonal

[m n]=size(U);
if(m~=2*n)
    error('cbrpack:oddSize','the input matrix must have an even number of rows');
end

v=false(n,1);
Q=eye(n);

%initialU=U;

firstPivot=nan;

Perm=1:n; %we permute the rows in each block to keep the already-factorized part in the first k-1 rows

for k=1:n
    SquaredRowNorms=sum(abs(U([k:n,n+k:2*n],k:end).^2),2);
    [pivotValue relativePivotRow]=max(SquaredRowNorms); %relative=relative to the "reduced" matrix

    if relativePivotRow>n-k+1
        p=relativePivotRow+2*(k-1)-n;
        U([p,n+p],:)=[U(n+p,:);-U(p,:)];
        v(Perm(p))=true;
    else
        p=relativePivotRow+(k-1);
    end
    
    %swaps row p and p+n into k and k+n
    U([p p+n k k+n],:)=U([k k+n p p+n],:);
    Perm([p k])=Perm([k p]);

%    [w,beta,rowNorm]=gallery('house',U(k,k:end).'); %H=I-beta*v*v'
%  we re-implement gallery/house since the function call took most of the
%  time here
    w=U(k,k:end).';
    if w(1)==0
        rowNorm=-norm(w);
    else
        rowNorm=-norm(w)*sign(w(1));
    end
    w(1)=w(1)-rowNorm;
    bbeta = -1/(rowNorm'*w(1));
    rowNorm=abs(rowNorm);
%
%
    U(:,k:end)=U(:,k:end)-bbeta*(U(:,k:end)*w)*w';
    U(k,k+1:end)=0;
    
    Q(k:end,:)=Q(k:end,:)-bbeta*w*(w'*Q(k:end,:));
    if isnan(firstPivot)
        firstPivot=rowNorm;
    end
    lastPivot=rowNorm;
    
%    tmp1=U*Q;
%    tmp2=rowSwap(initialU,v,'N');
%    norm(tmp1-tmp2([Perm,n+Perm],:))
end

%undoes the permutation
U(Perm,:)=U(1:n,:);
U(Perm+n,:)=U(n+1:2*n,:);

invcond=abs(lastPivot/firstPivot);
