function k=checkECAREInvariantSubspaceResidual(A,B,Q,R,S,U,stab)
% assess the quality of a stable invariant subspace
%
% k=checkCAREInvariantSubspaceResidual(A,B,Q,R,S,U,X,stab)
%
% A,B,Q,R,S determine the ECARE;
% U is the 2nxn stable invariant subspace
%
% stab can be 'stabilizing' (default) or 'antistabilizing'
%
% k.stability=how far the invariant subspace is from stability. 0=stable, large=bad
% k.pencilBackwardError=(normalized) 2-norm of the smallest perturbation to A-sE
%    that has the prescribed invariant subspace as exact invariant subspace
%    (see Byers, Benner, A structure-preserving method for generalized...)
% k.isGood=boolean value telling if the solution is "roughly ok"
% k.lagrangianityResidual = how fare we are from Lagrangianity
%

if not(exist('stab','var')) || isempty(stab)
    stab='stabilizing';
end

if(any(isnan(U)))
    k.residual=NaN;
    k.stability=NaN;
    k.riccatiResidual=NaN;
    return
end

[n m]=size(B);

if not(exist('S','var')) || isempty(S)
    S=zeros(size(B));
end

k.lagrangianityResidual=norm(U'*jay(2*n)*U);
[Uorth,useless]=qr(U,0);
Ue=[Uorth zeros(2*n,m);zeros(m,n) eye(m)]; %Ue is already orthogonal
[AA,EE]=evenPencil(A,B,Q,R,S);

[u s v]=svd([EE*Ue AA*Ue]);
s=diag(s);
k.pencilBackwardError=s(n+m+1)/s(1);

EAhat=diag(s(1:n+m))*v(:,1:n+m)';

Ehat=EAhat(:,1:n+m);Ahat=EAhat(:,n+m+1:end);
es=eig(Ehat,Ahat);

maxUnstable=max(real(es(isfinite(es))));
eigModulus=max(abs(es(isfinite(es))));

k.stability=maxUnstable/eigModulus;


k.isGood = (k.stability<1e-6) && k.pencilBackwardError<1e-10;
