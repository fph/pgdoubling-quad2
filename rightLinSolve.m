function [X invcond]=rightLinSolve(B,A,varargin)
% analogous of linsolve() for transposed linear systems
%
% [X invcond]=rightLinSolve(B,A,varargin)
%
% like linsolve(A,B,varargin), but solves the system X=B/A instead of X=A\B

[X invcond]=linsolve(A',B',varargin{:});
X=X';
