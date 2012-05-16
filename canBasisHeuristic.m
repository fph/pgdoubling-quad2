function p=canBasisHeuristic(U,varargin)
% find a good starting permutation p for the canonical basis
% namely, U(p(1:size(U,2)),:) shouldn't be too ill-conditioned
%
% p=canBasisHeuristic(U)
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling



[Q,R,p]=qr(U',0);
