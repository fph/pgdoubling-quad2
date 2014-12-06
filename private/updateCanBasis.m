function [X p invcond]=updateCanBasis(X,p,in,out)
% computes efficiently a new canonical basis from an old one 
%
% [newX newp invcond]=updateCanBasis(oldX,oldp,out,in)
%
% rows of I in the position(s) "out" are replaced with
% rows of X in the position(s) "in"
%
% warning: in and out here denote row indices in I and X, not in U
%
% you should call this in the optimization loop choosing in and out
% so that inv(X(in,out)) is small
%
% no error checking since this is meant to be called in a tight loop
% if you choose to call this directly, you're on your own

% the real update formula (obtained by a SMW trick) is
%newX=X+(xi+e_in)*inv(S)*(e_out-xj);
% OR X+(e_in+Xout)*inv(S)*(e_out'-e_in'*X)
% where e_in and e_out are part of the identity matrix
% however, this would lead to cancellation when S is very large
% therefore, we rearrange the computations as follow

n=size(X,2);
r=size(X,1);
i1=true(r,1);
i1(in)=false;
i2=in;

j1=true(n,1);
j1(out)=false;
j2=out;

invcond=norm(X(i2,j2));
X(i2,j2)=inv(X(i2,j2));
invcond=invcond*norm(X(i2,j2));

X(i2,j1)=-X(i2,j2)*X(i2,j1);
X(i1,j1)=X(i1,j1)+X(i1,j2)*X(i2,j1);
X(i1,j2)=X(i1,j2)*X(i2,j2);

p([n+in out])=p([out n+in]);

end
