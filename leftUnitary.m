function [H,Mnew] = leftUnitary(M,rows,cols)
% Compute a unitary matrix H such that Mnew = HM and the following holds:
% Mnew[rows==1, cols==1] is a zero matrix, Mnew[rows==0, cols==1] is square.

% (Selected columns of M are transformed to P[O; X], where P is a permutation and X square, by premultiplying M with a unitary H.)

% ROWS and COLS are 0/1 vectors, [length(rows),length(cols)] must be equal to size(M) -- not checked.

[H,Mnew] = rightUnitary(M',cols,rows);
H = H';
Mnew = Mnew';

end