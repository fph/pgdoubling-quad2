function X=vpaLyap(A,Q)
% vpa-capable version of lyap()
%
% X=vpaLyap(A,Q)
%
% works exactly as lyap(A,Q), but now A and Q can be variable precision
% arithmetic. Uses Kronecker product formulation, O(n^6): ok only for small
% data sets.
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

n=length(A);
X=-(kron(eye(n),A)+kron(A,eye(n)))\reshape(Q,n*n,1);
X=reshape(X,n,n);
