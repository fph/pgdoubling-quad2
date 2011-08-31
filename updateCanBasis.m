function [newX newp]=updateCanBasis(X,oldp,in,out)
% computes efficiently a new canonical basis from an old one 
%
% [newX newv]=updateCanBasis(oldX,oldv,out,in)
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

k=size(X,2);

S=X(in,out);
xj=X(in,:);
xi=X(:,out);

% the real update formula (obtained by a SMW trick) is
%newX=X+(xi+e_in)*inv(S)*(e_out-xj);
% OR X+(e_in+Xout)*inv(S)*(e_out'-e_in'*X)
% where e_in and e_out are part of the identity matrix
% however, this would lead to cancellation when S is very large
% therefore, we rearrange the computations as follow

%TODO: matgic:minv for this formula

xiTimesInvS=xi/S;

newX=X-xiTimesInvS*xj;

newX(in,:)=-S\xj;
newX(:,out)=xiTimesInvS;
newX(in,out)=inv(S);

newp=oldp;
[newp(k+in) newp(out)]=deal(newp(out),newp(k+in));

end
