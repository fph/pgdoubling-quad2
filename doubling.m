function [X,Y,U,V]=doubling(S,v,varargin)
% applies doubling to a symplectic pencil
%
% [X,Y,U,V]=doubling(X,v,varargin)
%
% X,v are the symBasis of a symplectic pencil
%
% options available:
% maxSteps: maximum number of step before returning (default:200)
% tolerance: tolerance at which to stop, if attainable (default: eps)
%   (the algorithm will stop anyway if the results aren't improving anymore)
% verbose: print convergence history and extra diagnostic information

o=matgic.Options(varargin{:});

maxSteps=o.get('maxSteps',inf);
tol=o.get('tolerance',eps);

verbose=logical(o.get('verbose',false));

steps=0;

convergenceHistory=[nan nan nan nan nan nan];

w=[]; %permutation guess
while(true)
    steps=steps+1;
    [S,v,w,swaps1,swaps2,nGH,nEF]=doublingStep(S,v,w);
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

n=length(S)/2;
first=1:n;second=n+1:2*n;

U=rowSwap([eye(n);-S(second,second);],v(second),'N');
X=U(second,:)/U(first,:);

V=rowSwap([-S(first,first);eye(n)],v(first),'T');
Y=V(first,:)/V(second,:);
