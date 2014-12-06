function wuv=wuvFromCanBasis(can)
% turns a canBasis into a wuv representation of a matrix
%
% wuv = wuvFromCanBasis(can)
%
% returns a structure wuv containing wuv.w, wuv.u, wuv.v such that M=w^T(I-uv^T)^{-1}, where
% M=canBasis2Matrix(X,p)
% wuv contains also wuv.ImvTu = I-v^Tu, computed in a more accurate way.

[n m]=size(can.X);
invp(can.p)=1:length(can.p);

%reorders columns so that the top block is identity + low rank
%Ystaying=find(invp(1:m)<=m);
Yleaving=find(invp(1:m)>m);
s=length(Yleaving);
Zentering=find(invp(m+1:end)<=m);

columnPermutation=invp(1:m);
columnPermutation(Yleaving)=invp(Zentering+m);
%I=eye(m);I(columnPermutation,:)
wuv.u=zeros(m,s);
wuv.u(Yleaving,:)=eye(length(Yleaving));

wcolumns=invp(m+1:end);
wuv.w=nan(n,m);
%populates w' with the rows from the identity
wuv.w(wcolumns<=m,:)=0;
wuv.w(wcolumns<=m,wcolumns(wcolumns<=m))=eye(sum(wcolumns<=m));
%populates w' with the rows from X
wuv.w(wcolumns>m,:)=can.X(wcolumns(wcolumns>m)-m,:);
wuv.w=wuv.w';

wuv.v=can.X(invp(Yleaving)-m,:);
wuv.v=wuv.v';

wuv.v=wuv.v(columnPermutation,:);
wuv.w=wuv.w(columnPermutation,:);

wuv.ImvTu=wuv.v'*wuv.u;
wuv.v=wuv.u-wuv.v;

