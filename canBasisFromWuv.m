function [can,invcond]=canBasisFromWuv(wuv)
% converts a wuv to a CanBasis
%
% [can,invcond]=canBasisFromWuv(wuv)
%
% a "wuv" is a representation of the form w'*inv(I-uv') for a short fat matrix.

[n m]=size(wuv.w);
r=size(wuv.u,2);

M=[eye(n)-wuv.u*wuv.v';wuv.w'];

canu=canBasisFromSubspace(wuv.u);

Y=wuv.u(canu.p(1:r),:);

wuv.v=productWithSparseL(canu,wuv.v','IT')';wuv.v=wuv.v([r+1:end 1:r],:);
wuv.w=productWithSparseL(canu,wuv.w','IT')';wuv.w=wuv.w([r+1:end 1:r],:);

L=productWithSparseL(canu,eye(n),'N');
L=L([r+1:end 1:r],:);
%v=L'\v;
%w=L'\w;
%MM=blkdiag(L,eye(m))*M*inv(L);

blk=[-Y*wuv.v';wuv.w'];blk(1:r,n-r+1:end)=eye(r)+blk(1:r,n-r+1:end); %TODO: get rid of that sum
Z2=null(blk(:,n-r+1:end)')';
Z1=-Z2*blk(:,1:n-r);
Z=[Z1 Z2];

%Z(:,1:n)=Z(:,1:n)*L;
Z(:,1:n)=Z(:,[n-r+1:n,1:n-r]);
Z(:,1:n)=productWithSparseL(canu,Z(:,1:n),'T');

can=canBasisFromSubspace(Z');
can.p=can.p([m+1:end 1:m]);

can.X=-can.X';
