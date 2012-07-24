function [w,u,v,ImvTu]=canBasis2wuv(X,p)
% turns a canBasis into a wuv representation of a matrix
%
% [w,u,v,ImvTu]=canBasis2wuv(X,p)
%
% returns w,u,v and I-v^Tu such that M=w^T(I-uv^T)^{-1}, where
% M=canBasis2Matrix(X,p)

[n m]=size(X);
invp(p)=1:length(p);

%reorders columns so that the top block is identity + low rank
%Ystaying=find(invp(1:m)<=m);
Yleaving=find(invp(1:m)>m);
s=length(Yleaving);
Zentering=find(invp(m+1:end)<=m);

columnPermutation=invp(1:m);
columnPermutation(Yleaving)=invp(Zentering+m);
%I=eye(m);I(columnPermutation,:)
u=zeros(m,s);
u(Yleaving,:)=eye(length(Yleaving));

wcolumns=invp(m+1:end);
w=nan(n,m);
%populates w' with the rows from the identity
w(wcolumns<=m,:)=0;
w(wcolumns<=m,wcolumns(wcolumns<=m))=eye(sum(wcolumns<=m));
%populates w' with the rows from X
w(wcolumns>m,:)=X(wcolumns(wcolumns>m)-m,:);
w=w';

v=X(invp(Yleaving)-m,:);
v=v';

v=v(columnPermutation,:);
w=w(columnPermutation,:);

ImvTu=v'*u;
v=u-v;

