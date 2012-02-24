function [A,B,Q,R,S]=hInfinityControlPencil(type,gamma,A,B1,B2,C1,C2,D11,D12,D21);
% returns one of the two pencils needed in H-infinity control
%
% [A,B,Q,R,S]=hInfinityControlPencil(type,gamma,A,B1,B2,C1,C2,D11,D12,D21);
%
% type == 'H' or 'J'

[n m1]=size(B1);
[n m2]=size(B2);
[p1 n]=size(C1);
[p2 n]=size(C2);


Q=zeros(n);

switch type
    case 'H'
        A=-A';
        B=[zeros(n,m1) zeros(n,m2) -C1'];
        S=[B1 B2 zeros(n,p1)];
        R=[gamma^2*eye(m1) zeros(m1,m2) D11'; zeros(m2,m1) zeros(m2) D12';D11 D12 eye(p1)];
    case 'J'
        A=-A;
        B=[zeros(n,p1) zeros(n,p2) -B1];
        S=[C1' C2' zeros(n,m1)];
        R=[gamma^2*eye(p1) zeros(p1,p2) D11; zeros(p2,p1) zeros(p2) D21;D11' D21' eye(m1)];
    otherwise
        error 'Wrong type'
end
