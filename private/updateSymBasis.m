function [X, newv, invcond]=updateSymBasis(X, oldv, inout)
% computes efficiently a new symplectic basis from an old one 
%
% [X, newv, invcond]=updateSymbBasis(X, oldv, inout)
%
% computes a new representation of the same subspace using a different v
% newv is the same as oldv, but with the bits at the indices in "inout" flipped
%
% it is effectively a PPT [Tsatsomeros, LAA '00] + some sign swaps
%
% no error checking since this is meant to be called in a tight loop
% if you choose to call this directly, you're on your own

w=X(:,inout);
invcond=norm(X(inout,inout));
s=inv(X(inout,inout));
invcond=invcond*norm(s);

X=X-w*s*w';
X(:,inout)=w*s;
X(inout,:)=X(:,inout)';
X(inout,inout)=-s;

toSwapSign= oldv(inout)==1;
X(:,inout(toSwapSign))=-X(:,inout(toSwapSign)); 
X(inout(toSwapSign),:)=-X(inout(toSwapSign),:); %TODO: can restructure to avoid this row?

newv=oldv;
newv(inout)=not(newv(inout));
