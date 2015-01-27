function[X, v, invcond] = updateQuadBasisIn(X, v, in)
% update a quadbasis by swapping a 0 to a 1
%
%  [X, newv, invcond] = updateQuadBasisIn(X, oldv, in)
%
% Computes a new quadBasis starting from an existing one. The new quadBasis
% needs to have newv(in)=true, while the old one needs to have oldv(in)=false.
% The rest of v is unchanged.
% in is an integer (index to ppt on).
%
% The behavior when oldv(in) is true is undefined.
% We do not check, and in general there is **no error checking**
% since this is meant to be called in a tight loop.
% If you choose to call this function directly, you're on your own.

% check updateQuadBasisOut.m for more detailed comments

first = v;
third = ~v; third(in) = false;

% Householder transformation to put zeros in X(in, third)
[w beta s] = gallery('house', [X(in,in) X(in, third)]');
scalar_product = X(third,in)*w(1) + X(third,third)*w(2:end,1);
X(third,third) = X(third,third) - scalar_product*beta*w(2:end,1)';
X(third,in) = X(third,in) - scalar_product*beta*w(1)';
X(in,third) = 0;
%X(in, in) = conj(s); this will be overwritten

% PPT
s = conj(s);
X(in, in) = -1/s;
X(in, first) = X(in, first) / s;
X(third,first) = X(third, first) - X(third,in)*X(in,first);
X(third, in) = X(third, in) / s;
X(first, in) = 0;
v(in) = true;

invcond = 1;


