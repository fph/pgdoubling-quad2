function testVanillaSDA

[A,G,Q]=carex(1);

[X,Y]=vanillaSDA(A,G,Q);

k=checkCAREInvariantSubspaceResidual(A,G,Q,X);
assert(k.isGood);
k=checkCAREInvariantSubspaceResidual(A,G,Q,inv(Y),[],'antistabilizing');
assert(k.isGood);
