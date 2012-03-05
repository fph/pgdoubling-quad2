function [X,v,invcond]=hamiltonianPencil2SymBasis(A,E,v)
%extracts symBases of Hamiltonian pencils and matrices
%
% [X,v,invcond]=hamiltonianPencil2SymBasis(A,E,v)
%
% use E=eye() or E=[] if you just have a Hamiltonian matrix

n=length(A);

if not(exist('E','var')) || isempty(E)
    E=eye(n);
end

if not(exist('v','var')) || isempty(v)
    v=[];
end

J=jay(n);
U=[E'; (A*J)'];

[X,v,invcond]=symplecticSubspace2SymBasis(U,v);
