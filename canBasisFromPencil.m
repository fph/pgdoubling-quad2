function [can,invcond,swaps]=canBasisFromPencil(A,E,varargin)
% Bounded canonical basis representation of a pencil
%
% [can,invcond,swaps]=pencil2CanBasis(A,E,options);
%
% input: A-xE pencil
%
% output:
% can: a canBasis of [E';A']
%
% Options: see subspace2CanBasis: 'permutation', 'threshold', 'maxswaps'
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

[can,invcond,swaps] = canBasisFromSubspace([E';A'], varargin{:});
can.origin = 'pencil';