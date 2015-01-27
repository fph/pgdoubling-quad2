function U=symplecticSubspaceFromQuadBasis(quad)
% unpacks the quadBasis representation of a quasidefinite symplectic subspace
%
% U=symplecticSubspaceFromQuadBasis(quad)
%
% given quad representing a symplectic subspace [I;X] with X\geq 0, return a matrix spanning it.
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

U = symplecticSubspaceFromSymBasis(symBasisFromQuadBasis(quad));