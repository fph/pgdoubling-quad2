function J=jay(n);
%returns J (the symplectic "swap everything" matrix)
%
% J=jay(n)
%
% returns a nxn matrix; n must be even

if(mod(n,2)~=0) error 'n must be even (try 2*n)';end
m=n/2;
J=zeros(2*m);
J(1:m,m+1:2*m)=eye(m);
J(m+1:2*m,1:m)=-eye(m);
