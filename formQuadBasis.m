function quad = formQuadBasis(C,A,B)
% Given the factors C, A, B that represent a quasidefinite matrix [-C*C A*; A BB*] form
% the quadBasis representation quad of the same matrix as quad.X = [Csquare 0; A Bsquare]
% and quad.v = [true .. true false .. false], where the now square factors
% Csquare and Bsquare are formed by padding possibly rectangular C and B
% with zeros. (The optimization algorithm works with quads of this form).

% B is tall and thin, C short and fat -- not checked.

[Br,Bc] = size(B);
[Cr,Cc] = size(C);

if (Br~=Bc), B = [B zeros(Br,Br-Bc)]; end % BorigBorig* = BnewBnew*
if (Cr~=Cc), C = [zeros(Cc-Cr,Cc); C]; end % Corig*Corig = Cnew*Cnew 
clear quad
quad.v = logical([ones(1,Cc) zeros(1,Br)]);
quad.X(quad.v,quad.v) = C;
quad.X(~quad.v,quad.v) = A;
quad.X(~quad.v,~quad.v) = B;
quad.X(quad.v,~quad.v) = 0;