function [X,Y,U,V]=solveCARE(A,G,Q,varargin)
% solves a continuous-time algebraic Riccati equation using doubling
%
% [X,Y,U,V]=solveCARE(A,G,Q,...)
%
% options:
%
% minSteps, maxSteps, tolerance, verbose: as in doubling.m
%
% type: either 'sda' or 'sign'
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

o=Options(varargin{:});

type=o.get('type','sda');

H=hamiltonian(A,G,Q);
scaling=norm(H);
[S,v]=hamiltonianPencil2SymBasis(H,scaling*eye(size(H)));
[S,v]=optimizeSymBasis(S,v);
switch type
    case 'sda'
        [AA,EE]=symBasis2HamiltonianPencil(S,v);
        
        %Cayley transform
        gamma=o.get('gamma',1.1*length(S)*max(max(abs(S)))); %this should ensure that gamma does not collide with some eigenvalues
        if not(gamma>0)
            error 'gamma must be positive'
        end
        
        [S,v]=symplecticPencil2SymBasis(AA+gamma*EE,AA-gamma*EE);
        [S,v]=optimizeSymBasis(S,v);
    case 'sign'
        if o.isSet('gamma')
            error 'specifying gamma makes sense only for sda, not for matrix sign'
        end

        %do nothing, S is already ok as it is
end

[S,v]=doubling(S,v,type,o);

switch type
    case 'sda'
        n=length(S)/2;
        first=1:n;second=n+1:2*n;
        
        Xpi=-S(second,second);vx=v(second);
        Ypi=-S(first,first);vy=v(first);
        
        U=rowSwap([eye(n);Xpi;],vx,'N');
        [X invcond1]=rightLinSolve(U(second,:),U(first,:));
        
        V=rowSwap([Ypi;eye(n)],vy,'T');
        [Y invcond2]=rightLinSolve(V(first,:),V(second,:));
    case 'sign'
        [A,E]=symBasis2HamiltonianPencil(S,v);
        
        %TODO: replace with sth else structure-preserving...
        n=size(S)/2;
        first=1:n;second=n+1:2*n;
        [u s v]=svd(A+E);U=v(:,second);
        [u s v]=svd(A-E);V=v(:,second);
        
        
        [X invcond1]=rightLinSolve(U(second,:),U(first,:));
        [Y invcond2]=rightLinSolve(V(first,:),V(second,:));
        
end
