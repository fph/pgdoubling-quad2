function [S,v]=doubling(S,v,varargin)
% applies doubling to a symplectic pencil
%
% [X,v]=doubling(X,v,varargin)
%
% X,v are the symBasis of a symplectic pencil
%
% options available:
% maxSteps: maximum number of step before returning (default:200)
% tolerance: tolerance at which to stop, if attainable (default: eps)
%   (the algorithm will stop anyway if the results aren't improving anymore)
% verbose: (logical) print convergence history and extra diagnostic information
% safer: (logical) uses a "safe" slower O(n^3) heuristic at every step

o=matgic.Options(varargin{:});

maxSteps=o.get('maxSteps',inf);
tol=o.get('tolerance',eps);

verbose=logical(o.get('verbose',false));
safer=logical(o.get('safer',false));

steps=0;

convergenceHistory=[nan nan nan nan nan nan];

w=[];
wguess=[]; %permutation guess
vguess=[];
while(true)
    steps=steps+1;
    if(safer)
        wguess=[];
        vguess=[];
    else
        wguess=w;
        vguess=v;
    end
    [S,v,w,swaps1,swaps2,nGH,nEF]=doublingStep(S,v,wguess,vguess);
    if(verbose)
        fprintf('Step %3d, subspace change measure %5.2e, pseudo-residual measure %5.2e, swaps 2*%d+%d\n',steps,nGH,nEF,swaps1,swaps2);
    end
    if nEF<=tol
        break;
    end
    convergenceHistory=[convergenceHistory(2:end) nEF];
    if nEF<1e-6 && nEF>=mean(convergenceHistory) %stagnation
        break
    end
    if(steps>maxSteps)
        break;
    end
end

