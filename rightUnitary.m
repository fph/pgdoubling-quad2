function [H,Mnew] = rightUnitary(M,rows,cols)
% Compute a unitary matrix H such that Mnew = MH and the following holds:
% Mnew[rows==1, cols==1] is a zero matrix, Mnew[rows==1, cols==0] is square.

% (Selected rows of M are transformed to [O X]P, where P is a permutation and X square, by postmultiplying M with a unitary H.)

% ROWS and COLS are 0/1 vectors, [length(rows),length(cols)] must be equal to size(M) -- not checked.

m = sum(rows); % number of rows for the zero block
n = sum(cols); % number of columns in the zero block

if length(cols)-n ~= m, error('Incompatible indices.'), end;

T = [M(rows==1, cols==0) M(rows==1, cols==1)];
[H,R] = qr(T');

Mnew(rows==1, cols==0) = R(1:m,1:m)'; % the triangular block from the R factor
Mnew(rows==1, cols==1) = 0;

T = [M(rows==0, cols==0) M(rows==0, cols==1)];
T = T*H;

Mnew(rows==0, cols==0) = T(:,1:m);
Mnew(rows==0, cols==1) = T(:,m+1:end);

% Fix the row and column order of H
T = H;
H(cols==0, cols==0) = T(1:m,1:m);
H(cols==1, cols==0) = T(m+1:end,1:m); 
H(cols==0, cols==1) = T(1:m,m+1:end); 
H(cols==1, cols==1) = T(m+1:end,m+1:end);

end

