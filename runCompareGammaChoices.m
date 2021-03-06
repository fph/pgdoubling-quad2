function runCompareGammaChoices(A,G,Q,gammas)
% solve a CARE with different gamma values for comparison, plot results
%
% runCompareGammaChoices(A,G,Q,gammas)
%
% gamma should be a vector
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

M=[];
for gamma=gammas
    H=hamiltonian(A,G,Q);
    s=svd([H+gamma*eye(size(H)) H-gamma*eye(size(H))]);
    
    [X,Y,U,V]=solveCARE(A,G,Q,'gamma',gamma,'safer',true,'maxSteps',50);
    k=checkCAREInvariantSubspaceResidual(A,G,Q,U);
    M=[M;k.residual k.riccatiResidual k.pencilBackwardError min(s)/max(s)];
end

loglog(gammas,M);
