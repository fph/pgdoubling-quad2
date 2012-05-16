function [X invcond]=rightLinSolve(B,A,varargin)
% analogous of linsolve() for transposed linear systems
%
% [X invcond]=rightLinSolve(B,A,varargin)
%
% like linsolve(A,B,varargin), but solves the system X=B/A instead of X=A\B
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

[X invcond]=linsolve(A',B',varargin{:});
X=X';
