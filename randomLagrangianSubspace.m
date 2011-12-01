function Q=randomLagrangianSubspace(twon)
% returns a random Lagrangian subspace
%
% Q=randomLagrangianSubspace(twon)
%
% returns an orthogonal basis of a random Lagrangian subspace, hopefully equidistributed

S=randomOrthosymplecticMatrix(twon);
Q=S(:,1:twon/2);
