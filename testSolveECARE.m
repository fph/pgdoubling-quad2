function testSolveECARE

[A,G,Q,X,parout,B,R]=carex(3);
[X,Y,U,V]=solveCARE(A,G,Q);

[Xe,Ye,Ue,Ve]=solveECARE(A,B,Q,R);

assertElementsAlmostEqual(subspace(Ue,abs(jay(length(U)))*U),0);
assertElementsAlmostEqual(subspace(Ve,abs(jay(length(U)))*V),0);

assertVectorsAlmostEqual(X,Xe);
assertVectorsAlmostEqual(Y,Ye);
