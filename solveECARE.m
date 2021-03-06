function [X,Y,U,V]=solveECARE(A,B,Q,R,S,varargin)
% solves a CARE given in 3x3 even block form, and computes invariant subspaces
%
% [X,Y,U,V]=solveECARE(A,B,Q,R,S,varargin)
%
% returns the first two blocks of the stable and unstable invariant subspaces
% of evenPencil(A,B,C,Q,R,S), U and V.
%
% Accepts options as in solveCARE
%
% X and Y are obtained by writing them in form [X;I] and [I;Y]
% respectively. In other words, X is the stabilizing solution of the Riccati equation
% with coefficients A-BR^(-1)C', -BR^(-1)B', Q-CR^(-1)C', and Y the one of
% its dual equation.
%
% Notice though that the order of the blocks in the ECARE and the CARE
% approach are different, so U and V obtained by this function and
% solveCARE differ by a factor abs(jay()), and we get [X;I] instead of
% [I;X]. This is a fault of the standard notations, not mine.
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

o=Options(varargin{:});

type=o.get('type','sda');

[n m]=size(B);

if not(exist('S','var'))
    S=[];
end

v=o.get('initialv',[]);

[AA,EE]=evenPencil(A,B,Q,R,S);

sym=symBasisFromEvenPencil(AA,EE,n,m,v);
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
    case 'sign'
        if o.isSet('gamma')
            error 'specifying gamma makes sense only for sda, not for matrix sign'
        end
        
        %do nothing, S is already ok as it is
    case 'signWithAnExtraUselessCayleyTransform'
        %to check if instabilities in sda are really due to the Cayley
        %forth and back
        
        [AA,EE]=hamiltonianPencilFromSymBasis(sym);
        
        %Cayley transform
        gamma=o.get('gamma',1.1*length(sym.X)*max(max(abs(sym.X)))); %this should ensure that gamma does not collide with some eigenvalues
        if not(gamma>0)
            error 'gamma must be positive'
        end
        sym=symBasisFromSymplecticPencil(AA+gamma*EE,AA-gamma*EE);
        gamma

        %now we undo the Cayley that we just did
        
        [AA,EE]=symplecticPencilFromSymBasis(sym);
        sym=symBasisFromHamiltonianPencil(AA+EE,AA-EE); %reverses the Cayley --- scaling shouldn't be needed?
        
        type='sign';
    otherwise
        error 'unknown type'
end

sym=doubling(sym,type,o);

switch type
    case 'sda'
        n=length(sym.X)/2;
        first=1:n;second=n+1:2*n;
        
        U=rowSwap([eye(n);-sym.X(second,second);],sym.v(second),'N');
        [X invcond1]=rightLinSolve(U(first,:),U(second,:));
        
        V=rowSwap([-sym.X(first,first);eye(n)],sym.v(first),'T');
        [Y invcond2]=rightLinSolve(V(second,:),V(first,:));
    case 'sign'
        [A,E]=symBasis2HamiltonianPencil(S,v);
        
        %to be replaced by sth else structure-preserving...
        n=length(sym.X)/2;
        first=1:n;second=n+1:2*n;
        [u s v]=svd(A+E);U=v(:,second);
        [u s v]=svd(A-E);V=v(:,second);
        
        
        [X invcond1]=rightLinSolve(U(first,:),U(second,:));
        [Y invcond2]=rightLinSolve(V(second,:),V(first,:));        
end
