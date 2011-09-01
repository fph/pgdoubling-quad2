function Y=rowSwap(U,v,transp)
% apply a symplectic row swap matrix
%
% Y=rowSwap(U,v,'N') or Y=rowSwap(U,v,'T')
%
% left-multiply the 2n x k matrix U by the "symplectic row swap matrix" given by
%
% J_v=[diag(1-v) diag(v);-diag(v) diag(1-v)]
%
% or its transpose (depending on the last argument 'N' or 'T').
% v is a logical (0-1) vector of length n;
% if v=false(n,1), then J_v=I; if v=true(n,1), then J_v=J.


v=logical(v);
n=length(v);
if not(2*n==size(U,1))
    error('cbrpack:wrongPermutationLength','the row swap vector has wrong length (actual: %d, expected (half of the symplectic matrix size): %g)',n,size(U,1)/2);
end

A=U(1:n,:);
B=U(n+1:end,:);

switch(transp)
    case 'N'
        [A(v,:),B(v,:)]=deal(B(v,:),-A(v,:));
    case 'T'
        [A(v,:),B(v,:)]=deal(-B(v,:),A(v,:));
    otherwise
        error('cbrpack:wrongTransposedness', 'parameter "transp" must be either ''N'' or ''T''');
end
Y=[A;B];
