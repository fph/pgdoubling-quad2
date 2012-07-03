function [v,w]=productWithSparseL(X,p,v,w)
% apply the permuted lower triangular matrix that turns a space into its canBasis
%
% [Lv,wL]=productWithSparseL(X,p,v,w)
%
% a canBasis can be seen as a matrix L=[I 0; -X I]*Pi' that turns U into
% [Y;0]. This function applies this matrix to the given matrices to get L*v and w*L.

[m n]=size(X);
v=v(p,:);
v(n+1:end,:)=v(n+1:end,:)-X*v(1:n,:);

w(:,1:n)=w(:,1:n)-w(:,n+1:end)*X;
w(:,p)=w;
