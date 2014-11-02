function X = matrixFromCanBasis(can)
% turns a canbasis into a matrix
% 
% M = matrixFromCanBasis(can);
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

% in fact, the matrix representation is a particular canBasis, so we just
% have to get it

if not(strcmp(can.origin,'matrix'))
    warning('PGDoubling:wrongOrigin','You are converting to a matrix a canbasis that did not originate from a matrix');
end

[m n] = size(can.X);
can = canBasisFromCanBasis(can,1:m+n);
X = can.X;
