function J=jay(n)
%returns J (the symplectic "swap everything" matrix)
%
% J=jay(n)
%
% returns a nxn matrix; n must be even
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

if(mod(n,2)~=0) error 'n must be even (try 2*n)';end
m=n/2;
J=[zeros(m) eye(m); -eye(m) zeros(m)];
