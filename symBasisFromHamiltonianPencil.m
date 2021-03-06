function [sym,invcond,swaps]=symBasisFromHamiltonianPencil(A,E,varargin)
%extracts symBases of Hamiltonian pencils and matrices
%
% [sym,invcond,swaps]=symBasisFromHamiltonianPencil(A,E,options)
%
% use E=eye() or E=[] if you just have a Hamiltonian matrix
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

o=Options(varargin{:});

n=length(A)/2;
if not(n==round(n))
    error 'This pencil does not have even size'
end

if not(exist('E','var')) || isempty(E)
    E=eye(2*n);
end

first = 1:n;
second = n+1:2*n;
U=[A(:,first) -E(:,first) -E(:,second) A(:,second)]';

[sym,invcond,swaps]=symBasisFromSymplecticSubspace(U,o);
sym.origin='hamiltonianPencil';
