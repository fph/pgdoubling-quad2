function Q=randomLagrangianSubspace(twon)
% returns a random Lagrangian subspace
%
% Q=randomLagrangianSubspace(twon)
%
% returns an orthogonal basis of a random Lagrangian subspace, hopefully equidistributed
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

S=randomOrthosymplecticMatrix(twon);
Q=S(:,1:twon/2);
