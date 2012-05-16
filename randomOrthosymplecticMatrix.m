function S=randomOrthosymplecticMatrix(twon)
% generates a random real orthogonal symplectic matrix
%
% S=randomOrthosymplecticMatrix(twon)
%
% generates a twon x twon matrix
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

if mod(twon,2)==1
    error('cbrpack:oddSize','symplectic matrices must have even size');
end

n=twon/2;
A=randn(n)+1i*randn(n);
[Q R]=qr(A);
Q=Q*diag(exp(2*pi*1i*rand(n,1)));
Re=real(Q);Im=imag(Q);
S=[Re Im; -Im Re];
