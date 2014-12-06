function can=canBasisTranspose(can);
% returns canBasis of M'
%
% canT=canBasisTranspose(can);
%
% returns a canBasis of M', where M=matrixFromCanBasis(can)

[m n]=size(can.X);


if not(strcmp(can.origin,'matrix'))
    warning('PGDoubling:wrongOrigin','You are trying to transpose a canBasis that did not originate from a matrix');
end

%can.origin='subspace';
%U=subspaceFromCanBasis(can);
%M=canBasis2Matrix(X,p);

can.X=-can.X';
can.p=can.p([n+1:end,1:n]);

% now can is a canBasis of the subspace ker U', where U was the original
% subspace/canBasis
%
% TODO: what we did above is in fact a dual - merge it with leftDual somehow?

swap([n+1:n+m,1:n])=1:n+m;
can.p=swap(can.p);

toChangeSign=[false(1,m) true(1,n)];
toChangeSign=toChangeSign(can.p);

can.X(:,toChangeSign(1:m))=-can.X(:,toChangeSign(1:m));
can.X(toChangeSign(m+1:end),:)=-can.X(toChangeSign(m+1:end),:);

