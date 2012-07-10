function [X,p,invcond]=wuv2CanBasis(w,u,v)
% converts a wuv to a CanBasis with given p
%
% [X,p,invcond]=wuv2SpecifiedCanBasis(w,u,v,p)
%
% a "wuv" is a representation of the form w'*inv(I-uv') for a short fat matrix.

[n m]=size(w);
r=size(u,2);

M=[eye(n)-u*v';w'];

[Xu,pu]=subspace2CanBasis(u);

Y=u(pu(1:r),:);

v=productWithSparseL(Xu,pu,v','IT')';v=v([r+1:end 1:r],:);
w=productWithSparseL(Xu,pu,w','IT')';w=w([r+1:end 1:r],:);

L=productWithSparseL(Xu,pu,eye(n),'N');
L=L([r+1:end 1:r],:);
%v=L'\v;
%w=L'\w;
%MM=blkdiag(L,eye(m))*M*inv(L);

blk=[-Y*v';w'];blk(1:r,n-r+1:end)=eye(r)+blk(1:r,n-r+1:end); %TODO: get rid of that sum
Z2=null(blk(:,n-r+1:end)')';
Z1=-Z2*blk(:,1:n-r);
Z=[Z1 Z2];

%Z(:,1:n)=Z(:,1:n)*L;
Z(:,1:n)=Z(:,[n-r+1:n,1:n-r]);
Z(:,1:n)=productWithSparseL(Xu,pu,Z(:,1:n),'T');

[X p]=subspace2CanBasis(Z');
p=p([m+1:end 1:m]);

X=-X';
