function [A,E]=pencilFromCanBasis(can)
% "unpack" a canonical basis returning a pencil
%
% [A,E]=canBasis2Subspace(can)
%
% inverse of pencil2CanBasis
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling


if not(strcmp(can.origin,'pencil'))
    warning('PGDoubling:wrongOrigin','You are converting to a pencil a canbasis that did not originate from a pencil');
end

can.origin = 'subspace';
EAt=subspaceFromCanBasis(can);

n=size(EAt,1)/2;

E=EAt(1:n,:)';
A=EAt(n+1:end,:)';