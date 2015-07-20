function X=symmetricMatrixFromSymBasis(symb)
% returns the matrix X such that U = [I;X]
%
% X = symBasis2SymBasis(symb)
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

symb = symBasisFromSymplecticSubspace(symplecticSubspaceFromSymBasis(symb),'swap',false(length(symb.v),1));
X = symb.X;
