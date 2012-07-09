function [v,w]=productWithSparseL(X,p,v,mode)
% apply the permuted lower triangular matrix that turns a space into its canBasis
%
% [Lv,wL]=productWithSparseL(X,p,v,mode)
%
% a canBasis can be seen as a matrix L=[I 0; -X I]*Pi' that turns U into
% [Y;0]. This function applies this matrix to the given matrix.
%
% if mode=='N', computes L*v
% if mode=='T', computes v*L
% if mode=='I', computes L\v
% if mode=='IT', computes v/L
%

[m n]=size(X);

switch mode
    case 'N'
        v=v(p,:);
        v(n+1:end,:)=v(n+1:end,:)-X*v(1:n,:);
    case 'T'
        v(:,1:n)=v(:,1:n)-v(:,n+1:end)*X;
        v(:,p)=v;
    case 'I'
        v(n+1:end,:)=v(n+1:end,:)+X*v(1:n,:);
        v(p,:)=v;
    case {'IT','TI'}
        v=v(:,p);
        v(:,1:n)=v(:,1:n)+v(:,n+1:end)*X;
    otherwise
        error 'Wrong mode'
end
