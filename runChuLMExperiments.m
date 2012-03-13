function results=runChuLMExperiments(experiments,varargin)

o=matgic.Options(varargin{:});

if not(exist('experiments','var'))
    experiments=1:33;
end

if ~o.isSet('maxSteps')
    o.set('maxSteps',200);
end

for i=experiments
    fprintf('[%d]',i);
    [A,G,Q]=ChuLM07Carex(i);
    [X,Y,U,V]=solveCARE(A,G,Q,o);
    k=checkCAREInvariantSubspaceResidual(A,G,Q,U);
    assertTrue(k.isGood);
    results{i}=k;
end
fprintf('\n');

for i=experiments disp([i results{i}.residual]), end
