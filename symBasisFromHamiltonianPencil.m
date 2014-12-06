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

n=length(A);

if not(exist('E','var')) || isempty(E)
    E=eye(n);
end

J=jay(n);
U=[E'; (A*J)'];

[sym,invcond,swaps]=symBasisFromSymplecticSubspace(U,o);
sym.origin='hamiltonianPencil';
