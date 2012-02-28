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

%initialU=U;

k1=0;k2=0; %{k1,k2}=length of the already factorized part in the {1st,2nd} block;

Q=eye(b);

firstPivot=nan;

while(true)
    
    rowNorms=sum(abs(U([k1+1:n,n+k1+1:2*n,2*n+k2+1:end],k1+k2+1:end).^2),2);
    
    [pivotValue relativePivotRow]=max(rowNorms); %relative=relative to the "reduced" matrix
    if isnan(firstPivot)
        firstPivot=pivotValue;
    end
    lastPivot=pivotValue;
    
    if relativePivotRow<=2*(n-k1)
        %pivotRow is in one of the first two blocks
        if relativePivotRow>n-k1
            absolutePivotRow=relativePivotRow+2*k1;
            conjugateRow=absolutePivotRow-n;
        else
            absolutePivotRow=relativePivotRow+k1;
            conjugateRow=absolutePivotRow+n;
        end
        
        rowToZeroOut=k1+1;

        %swaps row absolutePivot into k1+1
        U([rowToZeroOut,absolutePivotRow],:)=U([absolutePivotRow,rowToZeroOut],:);
        Pi([rowToZeroOut,absolutePivotRow],:)=Pi([absolutePivotRow,rowToZeroOut],:);
%        rowNorms([rowToZeroOut,absolutePivotRow],:)=rowNorms([absolutePivotRow,rowToZeroOut],:);

        %and row conjugateRow into n+k1+1
        U([n+rowToZeroOut,conjugateRow],:)=U([conjugateRow,n+rowToZeroOut],:);
        Pi([n+rowToZeroOut,conjugateRow],:)=Pi([conjugateRow,n+rowToZeroOut],:);
%        rowNorms([n+rowToZeroOut,conjugateRow],:)=rowNorms([conjugateRow,n+rowToZeroOut],:);


    else
        %pivotRow is in the third block
        absolutePivotRow=relativePivotRow+2*k1+k2;
        rowToZeroOut=2*n+k2+1;
        
        %swaps it into 2*n+k2+1
        U([rowToZeroOut,absolutePivotRow],:)=U([absolutePivotRow,rowToZeroOut],:);
        Pi([rowToZeroOut,absolutePivotRow],:)=Pi([absolutePivotRow,rowToZeroOut],:);
%        rowNorms([rowToZeroOut,absolutePivotRow],:)=rowNorms([absolutePivotRow,rowToZeroOut],:);


    end
    
    [v,beta]=gallery('house',U(rowToZeroOut,k1+k2+1:end).'); %H=I-beta*v*v'
    U(:,k1+k2+1:end)=U(:,k1+k2+1:end)-beta*U(:,k1+k2+1:end)*v*v';
    U(rowToZeroOut,k1+k2+2:end)=0;
    
    Q(k1+k2+1:end,:)=Q(k1+k2+1:end,:)-beta*v*v'*Q(k1+k2+1:end,:);

%    assertVectorsAlmostEqual(U*Q,initialU(Pi,:));
    
    if relativePivotRow<=2*(n-k1) %increment one of the two indices
        k1=k1+1;
    else
        k2=k2+1;
    end
    
    if k1+k2>=b-1
        break;
    end
end

invcond=lastPivot/firstPivot;
