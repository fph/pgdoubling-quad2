function results=runCarexExperiments(experiments)

if not(exist('experiments','var'))
    experiments=1:19;
end

for i=experiments
    fprintf('[%d]',i);
    [A,G,Q]=carex(i);
    if i==17
        needsSafer=true;
    else
        needsSafer=false;
    end
    [X,Y,U,V]=solveCARE(A,G,Q,'maxSteps',200,'safer',needsSafer);
    k=checkCAREInvariantSubspaceResidual(A,G,Q,U);
    assertTrue(k.isGood);
    results{i}=k;
end
fprintf('\n');

for i=1:length(results) disp(results{i}.residual), end
