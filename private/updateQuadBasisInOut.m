function[X, v, invcond] = updateQuadBasisInOut(X, v, in, out)
% update a quadbasis by swapping a 1 to a 0 and a 0 to 1
%
%  [X, newv, invcond] = updateQuadBasis1(X, oldv, in, out)
%
% Computes a new quadBasis starting from an existing one. The new quadBasis
% will have newv(out)=false and newv(in)=true, while the old one needs to
% have oldv(out)=true and oldv(in)=false.
% The rest of v is unchanged.
% out, in are integer indices.
%
% See updateQuadBasisOut for a more thorough explanation.
%
% The behavior when the assumptions aren't respected is undefined.
% We do not check, and in general there is **no error checking**
% since this is meant to be called in a tight loop.
% If you choose to call this function directly, you're on your own.

% we first make two Householder transformation, to introduce zeros in the blocks
% called "first" and "third" (for consistency with updateQuadBasisOut and
% updateQuadBasisIn).

first = v; first(out) = false;
third = ~v; third(in) = false;

% Householder transformation to put zeros in X(first, out)
[w beta s] = gallery('house', [X(out,out);X(first,out)]);
scalar_product = w(1)'*X(out,first)+w(2:end,1)'*X(first,first);
X(first, first) = X(first, first) - w(2:end,1)*beta*scalar_product;
X(out, first) = X(out, first) - w(1)*beta*scalar_product;
X(first,out) = 0;
X(out,out) = s; %we can omit this since it will be overwritten

% Householder transformation to put zeros in X(in, third)
[w beta s] = gallery('house', [X(in,in) X(in, third)]');
scalar_product = X(third,in)*w(1) + X(third,third)*w(2:end);
X(third,third) = X(third,third) - scalar_product*beta*w(2:end,1)';
X(third,in) = X(third,in) - scalar_product*beta*w(1)';
X(in,third) = 0;
X(in, in) = conj(s); %we can omit this since it will be overwritten

% the actual ppt update
gamma = X(in,in);
alpha = X(out,out);
beta = X(in,out);
Delta = hypot(alpha*gamma,beta);

X(third,first) = X(third,first) - X(third,[out,in])/Delta^2*[alpha*gamma^2 conj(beta); -alpha*beta*gamma alpha^2*gamma]*X([out,in],first);

% the two additional minuses wrt notes are there to correct signs after the
% in->out exchange
[X(out,first),X(in,first)] = deal(-(gamma^2*alpha*X(out,first) + conj(beta)*X(in,first)) / Delta^2,...
                    (-beta*X(out,first) + alpha*X(in,first)) / Delta);
[X(third,out) X(third,in)] = deal(-(gamma*X(third,out)-beta*X(third,in))/Delta,...
   (conj(beta)*X(third,out)+alpha^2*gamma*X(third,in))/Delta^2);

%...as is the omitted minus in X(out,in)
X(out,out) = -gamma/Delta;
X(in,in) = -alpha/Delta;
X(in,out) = 0;
X(out,in) = conj(beta)/Delta^2;

X(first,in) = 0;
X(out,third) = 0;

v(out) = false;
v(in) = true;

%TODO: invcond