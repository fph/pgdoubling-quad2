function [Pi,U,Q,invcond]=pirq(U)
%computes a "restricted" Pi*R*Q factorization of a subspace U
%
% [Pi,U,Q,invcond]=pirq(U)
%
% U has size 2n+m x n+m
%
% returns Pi, R, Q such that R*Q=U(Pi,:)
%
% Pi should "behave well" on the first two structured blocks
% Pi=blkdiag(Pi1,Pi1,Pi2), size(Pi1)==[n,n], size(Pi2)==[m,m]
%
% Here R([1:n 2*n+1:end]) is lower triangular, Q is orthogonal

[a b]=size(U);
n=a-b;
m=b-n;

Pi=(1:a)';

initialU=U;

k1=0;k2=0; %{k1,k2}=length of the already factorized part in the {1st,2nd} block;

Q=eye(b);

firstPivot=nan;

while(true)
%    Ubefore=round(U)
%    U([k1+1:n,n+k1+1:2*n,2*n+k2+1:end],k1+k2+1:end)
    rowNorms=sum(abs(U([k1+1:n,n+k1+1:2*n,2*n+k2+1:end],k1+k2+1:end).^2),2);
    
%    rowNorms
    
    [pivotValue relativePivotRow]=max(rowNorms); %relative=relative to the "reduced" matrix
    
    %k=row to zero out, p=pivot row
    if relativePivotRow<=2*(n-k1)
        %pivotRow is in one of the first two blocks
        if relativePivotRow>n-k1
            p=relativePivotRow+2*k1-n;
            U([p,n+p],:)=jay(2)*U([p,n+p],:);
            Pi([p n+p])=Pi([n+p p]);
        else
            p=relativePivotRow+k1;
        end
        
        k=k1+1;

        %swaps row p and p+n into k and k+n
        U([p p+n k k+n],:)=U([k k+n p p+n],:);
        Pi([p p+n k k+n],:)=Pi([k k+n p p+n],:);
    else
        %pivotRow is in the third block
        p=relativePivotRow+2*k1+k2;
        k=2*n+k2+1;
        
        %swaps it into 2*n+k2+1
        U([k,p],:)=U([p,k],:);
        Pi([k,p],:)=Pi([p,k],:);
%       rowNorms([rowToZeroOut,absolutePivotRow],:)=rowNorms([absolutePivotRow,rowToZeroOut],:);

%        'swap'
%        absolutePivotRow,rowToZeroOut

    end
    [v,beta,rowNorm]=gallery('house',U(k,k1+k2+1:end).'); %H=I-beta*v*v'
    U(:,k1+k2+1:end)=U(:,k1+k2+1:end)-beta*U(:,k1+k2+1:end)*v*v';
    U(k,k1+k2+2:end)=0;
    
    Q(k1+k2+1:end,:)=Q(k1+k2+1:end,:)-beta*v*v'*Q(k1+k2+1:end,:);
    if isnan(firstPivot)
        firstPivot=rowNorm;
    end
    lastPivot=rowNorm;

    assertVectorsAlmostEqual(U(Pi,:)'*[jay(2*n) zeros(2*n,m);zeros(m,2*n+m)]*U(Pi,:),zeros(n+m));
    assertVectorsAlmostEqual(abs(U*Q),abs(initialU(Pi,:)));
    
    if relativePivotRow<=2*(n-k1) %increment one of the two indices
        k1=k1+1;
    else
        k2=k2+1;
    end
    
    if k1+k2>=b %we had b before, but this would avoid the last swap
        break;
    end
end

invcond=abs(lastPivot/firstPivot);
