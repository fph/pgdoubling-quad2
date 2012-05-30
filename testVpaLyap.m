function testVpaLyap

reset(RandStream.getGlobalStream);

n=5;
A=randn(n);
A=A-1.1*max(eig(A))*eye(n); %creates a stable A
Q=randn(n);
Q=Q*Q';
A=vpa(A);
Q=vpa(Q);
X=vpaLyap(A,Q);
M=max(max(double(abs(A*X+X*A'+Q))));
assert(M<1e-30);
