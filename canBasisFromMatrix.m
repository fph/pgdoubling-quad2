function [can,invcond,swaps]=canBasisFromMatrix(X,varargin)
% Bounded canonical basis representation of a matrix
%
% [can,invcond]=canBasisFromMatrix(M,options);
%
% input: M matrix
%
% output:
% can: such that U(can.p(1:end),:) and [I;can.X] span the same subspace, where U=[I;M]
%
% Options: see subspace2CanBasis: 'permutation', 'threshold', 'maxswaps'
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

o=Options(varargin{:});

[m n]=size(X);
can.X = X;
can.p = 1:m+n;
can.origin = 'matrix';

if o.isSet('permutation')
    [can invcond]=canBasisFromCanBasis(can,o.get('permutation'));
    swaps=0;
else
    [can invcond swaps]=optimizeCanBasis(can,o.get('threshold',[]),o.get('maxswaps',[]));
end
