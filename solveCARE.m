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
sym=symBasisFromHamiltonianPencil(H,scaling*eye(size(H)));
sym=optimizeSymBasis(sym);
switch type
    case 'sda'
        [AA,EE]=hamiltonianPencilFromSymBasis(sym);
        
        %Cayley transform
        gamma=o.get('gamma',1.1*length(sym.X)*max(max(abs(sym.X)))); %this should ensure that gamma does not collide with some eigenvalues
        if not(gamma>0)
            error 'gamma must be positive'
        end
        
        sym=symBasisFromSymplecticPencil(AA+gamma*EE,AA-gamma*EE);
        sym=optimizeSymBasis(sym);
    case 'sign'
        if o.isSet('gamma')
            error 'specifying gamma makes sense only for sda, not for matrix sign'
        end

        %do nothing, S is already ok as it is
end

sym=doubling(sym,type,o);

switch type
    case 'sda'
        n=length(sym.X)/2;
        first=1:n;second=n+1:2*n;
        
        Xpi=-sym.X(second,second);vx=sym.v(second);
        Ypi=-sym.X(first,first);vy=sym.v(first);
        
        U=rowSwap([eye(n);Xpi;],vx,'N');
        [X invcond1]=rightLinSolve(U(second,:),U(first,:));
        
        V=rowSwap([Ypi;eye(n)],vy,'T');
        [Y invcond2]=rightLinSolve(V(first,:),V(second,:));
    case 'sign'
        [A,E]=hamiltonianPencilFromSymBasis(sym);
        
        %TODO: replace with sth else structure-preserving...
        n=size(sym.X)/2;
        first=1:n;second=n+1:2*n;
        [u, s, v]=svd(A+E);U=v(:,second);
        [u, s, v]=svd(A-E);V=v(:,second);
        
        
        [X, invcond1]=rightLinSolve(U(second,:),U(first,:));
        [Y, invcond2]=rightLinSolve(V(first,:),V(second,:));
        
end
