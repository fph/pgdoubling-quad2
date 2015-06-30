function [N,K,L,H] = inverseQuasidef(C,A,B)
% Compute the factors of a quasidefinite inverse of a quasidefinite matrix
% given in factored form.

% If M = [-C^*C A^*; A BB^*] with B and C square (can be of different size)
% then M^{-1} = [-NN^* K; K^* L^*L];
% H is a unitary matrix whose blocks are needed for quadBasisFromQuadBasis,
% Case 3.

[m,n] = size(A); % C is nxn, B is mxm

[Q,R] = qr([B';A']);
R = R(1:m,:);

M = [Q(1:m,:); C*Q(m+1:end,:)];
[H,Mnew] = leftUnitary(M,[true(1,m) false(1,n)],[false(1,m) true(1,n)]);
H=H';

[N, invcondN] = rightLinSolve(Q(m+1:end,m+1:end),Mnew(m+1:end,m+1:end)); % N = Q_{22}Mnew_{22}^{-1} = BA^{-1}
[L, invcondL] = rightLinSolve(Mnew(1:m,1:m),R'); % L = Mnew_{11}R^{-*}

X = Q(m+1:end,1:m) - N*Mnew(m+1:end,1:m);
K = rightLinSolve(X,R');

if invcondN < sqrt(eps)
    warning('inverseQuasidef:illConditionedMatrix', 'The matrix Mnew_{22} has conditioning >1/sqrt(eps).');
end

if invcondL < sqrt(eps)
    warning('inverseQuasidef:illConditionedMatrix', 'The matrix R has conditioning >1/sqrt(eps).');
end

end

