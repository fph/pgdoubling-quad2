function M = matrixFromQuadBasis(quad)
% converts a quadBasis to its generating matrix

%fills up top B block

M = quad.X;
M(quad.v,~quad.v)=M(~quad.v,quad.v)';
M(quad.v,quad.v) = -M(quad.v,quad.v)'*M(quad.v,quad.v);
M(~quad.v,~quad.v) = M(~quad.v,~quad.v)*M(~quad.v,~quad.v)';