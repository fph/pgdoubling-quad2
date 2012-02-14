function results=runCarexExperiments(experiments)

if not(exist('experiments','var'))
    experiments=1:19;
end

for i=experiments
    [A,G,Q]=carex(i);
    [X,Y,U,V]=solveCARE(A,G,Q);
    k=checkCAREInvariantSubspaceResidual(A,G,Q,U);
    assertTrue(k.isGood);
    results{i}=k;
end
