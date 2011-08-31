function testSubspace2CanBasis

U=[1 1;1 1; 0 4];

[X,p]=subspace2CanBasis(U);
U2=canBasis2Subspace(X,p);

assertVectorsAlmostEqual(subspace(U,U2),0);

U=[1 1; 1 1; 1 1+1e-10];

assertWarningThrown(@() canBasisHeuristic(U),'cbrpack:illConditionedSubspace');

