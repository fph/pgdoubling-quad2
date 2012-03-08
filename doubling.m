function [S,v]=doubling(S,v,type,varargin)
% applies doubling to a symplectic or inverse-free sign to a Hamiltonian pencil
%
%  [X,v]=doubling(X,v,type,varargin)
%
% if type=='sign', (X,v) is a symBasis of a Hamiltonian pencil, and the
% algorithm performs repeated H <- H+inv(H)
%
% if type=='sda', (X,v) is a symBasis of a symplectic pencil, and the
% algorithm performs repeates S <- S*S
%
%
% options available:
% maxSteps: maximum number of step before returning (default:200)
% minSteps: minimum number of steps
% tolerance: tolerance at which to stop, if attainable (default: eps)
%   (the algorithm will stop anyway if the results aren't improving anymore)
% verbose: (logical) print convergence history and extra diagnostic information
% safer: (logical) uses a "safe" but slower QRP factorization at each step

o=matgic.Options(varargin{:});

maxSteps=o.get('maxSteps',inf);
minSteps=o.get('minSteps',0);
tol=o.get('tolerance',eps);

verbose=logical(o.get('verbose',false));
safer=logical(o.get('safer',false));

switch type
    case 'sda'
        f=@doublingStep;
    case 'sign'
        f=@inverseFreeSignStep;
    otherwise
        error 'invalid doubling type selected'
end

steps=0;

convergenceHistory=[nan nan nan nan nan nan];

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
    [S,v,w,swaps1,swaps2,res,res2]=f(S,v,wguess,vguess);
    if(verbose)
        fprintf('Step %3d, residual 1 %5.2e, residual 2 %5.2e, swaps 2*%d+%d\n',steps,res,res2,swaps1,swaps2);
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
