function testSolveCare

%check all parts of the solution on "easy" examples

for i=1:9
    [A,G,Q]=carex(i);
    [X,Y,U,V]=solveCARE(A,G,Q);
    k=checkCAREInvariantSubspaceResidual(A,G,Q,X);
    assertTrue(k.isGood);
    k=checkCAREInvariantSubspaceResidual(A',-Q,-G,Y);
    assertTrue(k.isGood);
    k=checkCAREInvariantSubspaceResidual(A,G,Q,U);
    assertTrue(k.isGood);
    k=checkCAREInvariantSubspaceResidual(A,G,Q,V,[],'antistabilizing');
    assertTrue(k.isGood);    
end

%check only what was meant to be solved in carex for the others

for i=10:19
    [A,G,Q]=carex(i);
    [X,Y,U,V]=solveCARE(A,G,Q);
%    k=checkCAREInvariantSubspaceResidual(A,G,Q,X);
%    assertTrue(k.isGood);
%    k=checkCAREInvariantSubspaceResidual(A',-Q,-G,Y);
%    assertTrue(k.isGood);
    k=checkCAREInvariantSubspaceResidual(A,G,Q,U);
    assertTrue(k.isGood);
%    k=checkCAREInvariantSubspaceResidual(A,G,Q,V,[],'antistabilizing');
%    assertTrue(k.isGood);    
end
