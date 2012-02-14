function [X,Y,U,V]=solveCARE(A,G,Q,varargin)
% solves a continuous-time algebraic Riccati equation using permuted doubling
%
% [X,Y,U,V]=solveCARE(A,G,Q,...)
%
% options:
% maxSteps: maximum number of step before returning (default:200)
% tolerance: tolerance at which to stop, if attainable (default: eps)
%   (the algorithm will stop anyway if the results aren't improving anymore)
%
% gamma: constant to use for the Cayley transform
% verbose: be verbose

o=matgic.Options(varargin{:});

H=hamiltonian(A,G,Q);
gamma=o.get('gamma',norm(H));
if not(gamma>0)
    error 'gamma must be positive'
end
maxSteps=o.get('maxSteps',inf);
tol=o.get('tolerance',eps);

verbose=logical(o.get('verbose',false));

[S,v]=symplecticPencil2SymBasis(H+gamma*eye(size(H)),H-gamma*eye(size(H)));

steps=0;

threshold1=2;
options2=matgic.Options();
%TODO: set options

w=[]; %permutation guess
while(true)
    steps=steps+1;
    options2.set('initialRowSwap',v);
    [S,v,w,swaps1,swaps2,nn]=doublingStep(S,v,w,threshold1,options2);
    if(verbose)
        fprintf('Step %3d, subspace change measure %5.2e, swaps 2*%d+%d\n',steps,nn,swaps1,swaps2);
    end
    if nn<=tol
        break;
    end
    if(steps>maxSteps)
        break;
    end
end

n=length(A);
first=1:n;second=n+1:2*n;

U=rowSwap([eye(n);-S(second,second);],v(second),'N');
X=U(second,:)/U(first,:);

V=rowSwap([-S(first,first);eye(n)],v(first),'T');
Y=V(first,:)/V(second,:);
