function v=productWithSparseL(can,v,mode)
% apply the permuted lower triangular matrix that turns a space into its canBasis
%
% w=productWithSparseL(can,v,mode)
%
% a canBasis can be seen as a matrix L=[I 0; -X I]*Pi' that turns U into
% [Y;0]. This function applies this matrix to the given matrix.
%
% if mode=='N', computes L*v
% if mode=='T', computes v*L
% if mode=='I', computes L\v
% if mode=='IT', computes v/L
%

[m n]=size(can.X);

switch mode
    case 'N'
        v=v(can.p,:);
        v(n+1:end,:)=v(n+1:end,:)-can.X*v(1:n,:);
    case 'T'
        v(:,1:n)=v(:,1:n)-v(:,n+1:end)*can.X;
        v(:,can.p)=v;
    case 'I'
        v(n+1:end,:)=v(n+1:end,:)+can.X*v(1:n,:);
        v(can.p,:)=v;
    case {'IT','TI'}
        v=v(:,can.p);
        v(:,1:n)=v(:,1:n)+v(:,n+1:end)*can.X;
    otherwise
        error 'Wrong mode'
end
