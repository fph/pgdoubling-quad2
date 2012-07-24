function [X,p]=canBasisTranspose(X,p);
% returns canBasis of M'
%
% [Xt,pt]=canBasisTranspose(X,p);
%
% returns a canBasis of M', where M=canBasis2Matrix(X,p)

[m n]=size(X);

U=canBasis2Subspace(X,p);
M=canBasis2Matrix(X,p);

X=-X';
p=p([n+1:end,1:n]);

% now X,p is a canBasis of the subspace ker U', where U was the original
% subspace/canBasis
%
% TODO: what we did above is in fact a dual - merge it with leftDual somehow?

swap([n+1:n+m,1:n])=1:n+m;
p=swap(p);

toChangeSign=[false(1,m) true(1,n)];
toChangeSign=toChangeSign(p);

X(:,toChangeSign(1:m))=-X(:,toChangeSign(1:m));
X(toChangeSign(m+1:end),:)=-X(toChangeSign(m+1:end),:);

