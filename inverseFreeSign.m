function sym=inverseFreeSign(sym,varargin)
% applies the inverse-free sign method to a Hamiltonian pencil
%
% sym=doubling(sym,varargin)
%
% sym is the symBasis of a Hamiltonian pencil
%
% options available:
% maxSteps: maximum number of step before returning (default:200)
% minSteps: minimum number of steps, independently of the stopping
% criterion
% tolerance: tolerance at which to stop, if attainable (default: eps)
%   (the algorithm will stop anyway if the results aren't improving anymore)
% verbose: (logical) print convergence history and extra diagnostic information
% safer: (logical) uses a "safe" slower O(n^3) heuristic at every step
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

o=Options(varargin{:});

maxSteps=o.get('maxSteps',inf);
minSteps=o.get('minSteps',0);
tol=o.get('tolerance',eps*length(S));

verbose=logical(o.get('verbose',false));
safer=logical(o.get('safer',false));

steps=0;

convergenceHistory=[nan nan nan nan nan nan];

if not(strcmp(sym.origin,'hamiltonianPencil'))
    warning('PGDoubling:wrongOrigin','You are converting to a Hamiltonian pencil a symBasis that did not come from one.');
end

w=[];
while(true)
    steps=steps+1;
    if(safer)
        wguess=[];
        vguess=[];
    else
        wguess=w;
        vguess=v;
    end
    [sym,w,swaps1,swaps2,res]=inverseFreeSignStep(sym,wguess,vguess);
    if(verbose)
        fprintf('Step %3d, pseudo-residual measure %5.2e, swaps 2*%d+%d\n',steps,res,swaps1,swaps2);
    end
    if steps<minSteps
        continue;
    end
    if res<=tol
        break;
    end
    convergenceHistory=[convergenceHistory(2:end) res];
    if res<1e-6 && res>=mean(convergenceHistory) %stagnation
        break
    end
    if(steps>=maxSteps)
        break;
    end
end
