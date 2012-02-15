function k=checkCAREInvariantSubspaceResidual(A,G,Q,U,X,stab)
% assess the quality of a stable invariant subspace
%
% k=checkCAREInvariantSubspaceResidual(A,G,Q,U,X,stab)
%
% A,G,Q determine the CARE;
% U is the 2nxn stable invariant subspace and X is the nxn Riccati solution
% One of them can be omitted; in this case it is computed from the other
% (U=[I;X] or X=U2/U1, respectively)
%
% stab can be 'stabilizing' (default) or 'antistabilizing'
%
% k.residual=subspace relative residual (as in ChuLM07, figure 4)
% k.stability=how far the invariant subspace is from stability. 0=stable, large=bad
% k.riccatiResidual=residual of X in the Riccati equation
% k.pencilBackwardError=(normalized) 2-norm of the smallest perturbation to A-sE
%    that has the prescribed invariant subspace as exact invariant subspace
%    (see Byers, Benner, A structure-preserving method for generalized...)
% k.isGood=boolean value telling if the solution is "roughly ok"
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

n=length(A);

if not(exist('X','var')) || isempty(X)
    %we have to reconstruct either U or X
    if(length(U)==n)
        %reconstruct U from X
        X=U;
        U=[eye(size(U));U];
    else
        %reconstruct X from U
        assertEqual(size(U),[2*n n]);
        [X invcond]=linsolve(U(1:n,:)',U(n+1:end,:)'); %this form does not throw warnings if the matrix to invert is ill-conditioned
        X=X';
    end
end

k.riccatiResidual=norm(A'*X+X*A+Q-X*G*X)/(norm(A'*X)+norm(X*A)+norm(Q)+norm(X*G*X));

k.lagrangianityResidual=norm(U'*jay(2*n)*U);

H=hamiltonian(A,G,Q);
[Uorth,useless]=qr(U,0);

[u s v]=svd([U H*U]);
s=diag(s);
k.pencilBackwardError=s(n+1)/s(1);

k.residual=norm(H*Uorth-Uorth*(Uorth'*H*Uorth))/norm(H); %they use norm 2
L=eig(Uorth'*H*Uorth);
switch stab
    case 'stabilizing'
        %do nothing
    case 'antistabilizing'
        L=-L;
    otherwise
        error 'wrong stability type requested'
end
emax=max(real(L));
k.stability=max(0,emax/norm(H));
k.stability=abs(k.stability);

k.isGood = (k.stability<1e-6) && k.residual<1e-10;
