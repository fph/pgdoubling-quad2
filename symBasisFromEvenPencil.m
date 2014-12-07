function [sym,invcond]=symBasisFromEvenPencil(AA,EE,n,m,v)
%converts an even pencil to a symBasis
%
% [sym,invcond]=symBasisFromEvenPencil(AA,EE,n,m,v)
%
% since we wish to allow more general pencils AA-sEE, we give n and m
% explicitly
%
% Warning: kinda misnamed, does not work with every even pencil --- I need
% to write down the theory
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

firstSecond=1:2*n;third=2*n+1:2*n+m;
U=[EE(:,firstSecond)'; (AA(:,firstSecond)*jay(2*n))'];

if not(exist('v','var')) || isempty(v)
    [v invcond1]=extendedSymBasisHeuristicPaper([U; AA(:,third)']);
end

U=rowSwap(U,v,'N');


[Xext invcond]=rightLinSolve(U(2*n+1:4*n,:),[U(1:2*n,:); AA(third,:)]);

sym.X=Xext(:,firstSecond);

%final check
if norm(sym.X-sym.X','fro')/norm(sym.X,'fro') > sqrt(eps)
    warning('cbrpack:notSymplectic','the resulting matrix is numerically very far from Hermitian --- was your starting subspace symplectic?');
end

%symmetrize
sym.X=(sym.X+sym.X')/2;
sym.v = v;
sym.origin='hamiltonianPencil'; %an even pencil "is" a Hamiltonian one, in this context

if(invcond<sqrt(eps(class(U))))
    warning('cbrpack:illConditionedMatrix', 'symplecticSubspace2SymBasis: the matrix I am inverting has conditioning >1/sqrt(eps). This may be due to an ill-conditioned subspace or to a bad initial guess --- consider using the initial value heuristic instead');
end
