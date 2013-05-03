function [X,p,invcond,swaps]=pencil2CanBasis(A,E,varargin)
% Bounded canonical basis representation of a pencil
%
% [X,p,invcond,swaps]=pencil2CanBasis(A,E,options);
%
% input: A-xE pencil
%
% output:
% X,p: a canBasis of [E';A']
%
% Options: see subspace2CanBasis: 'permutation', 'threshold', 'maxswaps'
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

[X,p,invcond,swaps]=subspace2CanBasis([E';A'],varargin{:});