function [X,Y]=vanillaSDA(A,G,Q,varargin)
%solves a CARE with the original SDA
%
% [X,Y]=vanillaSDA(A,G,Q)
%

o=matgic.Options(varargin{:});

verbose=logical(o.get('verbose',false));

maxIters=o.get('maxSteps',100);

H=hamiltonian(A,G,Q);
gamma=o.get('gamma',norm(H,'fro'));

n=length(A);

EE=H-gamma*eye(2*n);
AA=H+gamma*eye(2*n);

first=1:n;second=n+1:2*n;

EFGH=[EE(:,first) AA(:,second)]\[AA(:,first) EE(:,second)];

E=EFGH(first,first);
F=EFGH(second,second);
G=-EFGH(first,second);
H=-EFGH(second,first);

oldres=inf;
%S=[-G E; F -H];
for k=1:maxIters
%    Snew=symSDAStep(S);Gnew=-S(first,first);Hnew=-S(second,second);Enew=S(first,second);S=Snew;
    [L,U]=lu(eye(n) -G*H);
    Estar=E/U/L;
    Gnew=G+Estar*G*E';
    Hnew=H+(E'/L'/U')*H*E;
    Enew=Estar*E;
    %newres=norm(Gnew-G,'fro');
    newres=norm(Enew,'fro');
    cd=cond(eye(n)-G*H);
    if verbose
        fprintf('step %3d residual measures: %e %e -- cond %e\n',k,norm(Enew,'fro'),norm(Gnew-G,'fro'),cd)
    end
    if(newres>=oldres && newres<1.e-1)
        break;
    end
    E=Enew;G=Gnew;H=Hnew;
    oldres=newres;
end

X=H;
Y=G;
