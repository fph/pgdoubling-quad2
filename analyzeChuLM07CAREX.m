function analyzeChuLM07CAREX

% Check which experiments from
% D. Chu, X. Liu, and V. Mehrmann. A numerical method for computing the Hamiltonian
% Schur form, 2007
% can be used to define a quasidefinite matrix X which we can then use in
% the optimization algorithm.

% This is the "extended version" of the experiments in
% Benner, Laub, Mehrmann. A Collection of Benchmar Examples for the
% Numerical Solution of Algebraic Riccati Equations I: Continuous-Time Case

experiments = 1:33;
canUseBoth = [];
canUseFactors = [];

for i=experiments

    fprintf('[%d]: ',i);
    [A,G,Q,X, parout, B, R, C, Q0]=ChuLM07Carex(i);
    % these factors define X = [-Q A*; A G] = [-C*Q0C A*; A BR^{-1}B*]
       
    % we can define the factors of our Gram matrices BfBf* and Cf*Cf either
    % from the Cholesky decompositions of Q and G or from using the factors
    % B and C which define the factorizations of Q and G as above.
    [cholG,pG] = chol(G,'lower'); % G = BfBf* for us so we take Bf to be the lower triangular Cholesky factor, if we can
    [cholQ,pQ] = chol(Q); % Q = Cf*Cf for us so we take Cf to be the upper triangular Cholesky factor, if we can
    [cholQ0,pQ0] = chol(Q0); % Q = (cholQ0 C)*(cholQ0 C) =: Cf*Cf, if Q0 pos def
    [cholR,pR] = chol(R); % G=BR^{-1}B* = (B cholR^{-1})(B chol(R)^{-1})* =: BfBf*, if R pos def

    if pG, fprintf(1,'Cholesky G fails, min(eig(G)) = %5.3e ',min(eig(G))), end
    if pQ, fprintf(1,'Cholesky Q fails. min(eig(Q)) = %5.3e ',min(eig(Q))'), end
    if pQ0, fprintf(1,'Cholesky Q0 fails. min(eig(Q0)) = %5.3e ',min(eig(Q0))'), end
    if pR, fprintf(1,'Cholesky R fails. min(eig(R)) = %5.3e ',min(eig(R))'), end
    if [pG pQ pQ0 pR] == [0 0 0 0]
        canUseBoth = [canUseBoth i];
        fprintf(1,'Can use Choleskys of Q and G, or use the factors C and B.')
    end
    if [pQ0 pR] == [0 0]
        canUseFactors = [canUseFactors i];
    end
    fprintf(1,'\n')

end

fprintf(1,'\nFrom the above, we see that the experiments [2], [3], [4], [17] and [18] are problematic ')
fprintf(1,'for the reduction; \nby inspection, we see that the matrix Q0 in [2] is singular positive semidefinite and \n')
fprintf(1,'we can use cholQ0 = [3 2; 0 0] for this example so we add it to canUseFactors.\n')
fprintf(1,'Q0 in other experiments is indefinite so we cannot use them:\n')

[A,G,Q,X,parout,B,R,C,Q0]=ChuLM07Carex(3); fprintf(1,'[3], lmin(Q0) = %e\n',min(eig(Q0)))
[A,G,Q,X,parout,B,R,C,Q0]=ChuLM07Carex(4); fprintf(1,'[4], lmin(Q0) = %e\n',min(eig(Q0)))
[A,G,Q,X,parout,B,R,C,Q0]=ChuLM07Carex(17); fprintf(1,'[17], lmin(Q0) = %e\n',min(eig(Q0)))
[A,G,Q,X,parout,B,R,C,Q0]=ChuLM07Carex(18); fprintf(1,'[18], lmin(Q0) = %e\n',min(eig(Q0)))

canUseBoth
canUseFactors = sort([canUseFactors 2])