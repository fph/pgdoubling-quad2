function [AA,EE]=evenPencil(A,B,Q,R,S)
%constructs the even matrix pencil associated to a control problem
%
% [AA,EE]=evenPencil(A,B,Q,R,S)
%
%    [0  A  B]    [ 0  I 0]
% AA=[A' Q  S] EE=[-I  0 0]
%    [B' S' R]    [ 0  0 0]
%
% S may be omited (it is often zero)
% A,B,Q,R may be taken as the output of carex

if not(exist('S','var')) || isempty(S)
    S=zeros(size(B));
end

[n m]=size(B);
AA=[zeros(n) A B; A' Q S; B' S' R]; 
EE=zeros(size(AA));n=length(A); EE(1:n,n+1:2*n)=eye(n);EE(n+1:2*n,1:n)=-eye(n);
