function [X,p]=matrix2CanBasis(X)
% compute a representation E^-1A of a matrix
%
% [X,p]=matrix2CanBasis(X)

%internally, it is implemented through transposition

X=X';
n=size(X,2);
r=size(X,1);
p=1:n+r;
