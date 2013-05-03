function [X,p,invcond,swaps]=matrix2CanBasis(X,varargin)
% Bounded canonical basis representation of a matrix
%
% [X,p,invcond]=subspace2CanBasis(M,options);
%
% input: M matrix
%
% output:
% X,p: such that U(p(1:end),:) and [I;X] span the same subspace, where U=[I;M]
%
% Options: see subspace2CanBasis: 'permutation', 'threshold', 'maxswaps'
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

o=Options(varargin{:});

[m n]=size(X);
p=1:m+n;

if o.isSet('permutation')
    [X p invcond]=canBasis2CanBasis(X,p,o.get('permutation'));
    swaps=0;
else
    [X p invcond swaps]=optimizeCanBasis(X,p,o.get('threshold',[]),o.get('maxswaps',[]));
end
