function[X, v, invcond] = updateQuadBasisOut(X, v, out)
% update a quadbasis by swapping a 1 to a 0
%
%  [X, newv, invcond] = updateQuadBasis1(X, oldv, out)
%
% Computes a new quadBasis starting from an existing one. The new quadBasis
% needs to have newv(out)=false, while the old one needs to have oldv(out)=true.
% The rest of v is unchanged.
% out is an integer (index to ppt on).
%
% The behavior when oldv(out) is false is undefined.
% We do not check, and in general there is **no error checking**
% since this is meant to be called in a tight loop.
% If you choose to call this function directly, you're on your own.

first = v; first(out) = false;
third = ~v;

% Householder transformation to put zeros in X(first, out)
% this transformation should leave the subspace represented unaltered.
[w beta s] = gallery('house', [X(out,out);X(first,out)]);
scalar_product = w(1)'*X(out,first)+w(2:end,1)'*X(first,first);
X(first, first) = X(first, first) - w(2:end,1)*beta*scalar_product;
X(out, first) = X(out, first) - w(1)*beta*scalar_product;
X(first,out) = 0;
%X(out,out) = s; we omit this since it will be overwritten almost
%immediately

% real PPT update
% note that we are pivoting *back* 1->0 in a symplectic PPT, so there is
% an additional minus sign in row and column out
X(out, out) = - 1 / s;
X(out, first) = -X(out, first) / s;
X(third, first) = X(third, first) + X(third, out)*X(out, first); %the minus and the division are already in the previous line
X(third, out) = -X(third, out) / s;
X(out, third) = 0;
v(out) = false;

invcond = 1;