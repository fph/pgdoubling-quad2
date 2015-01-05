function quad = symBasisFromQuadBasis(quad)
% converts a quadBasis to a symBasis with the same v
%
% sym = symBasisFromQuadBasis(quad)
%
% whenever matrixFromSymBasis(sym) is positive definite, we can factor
% sym.X as [-A'*A B'; B C*C']. A /quadBasis/ is a structure containing 
% quad.v=sym.v and quad.X=[A arb; B C], where
% the block "arb" is arbitrary (never used, it's just there to fill in the matrix).

quad.X(quad.v,~quad.v)=quad.X(~quad.v,quad.v)';
quad.X(quad.v,quad.v) = -quad.X(quad.v,quad.v)'*quad.X(quad.v,quad.v);
quad.X(~quad.v,~quad.v) = quad.X(~quad.v,~quad.v)*quad.X(~quad.v,~quad.v)';
