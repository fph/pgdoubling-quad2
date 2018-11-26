function reduceChuLM07CAREX(threshold,varargin)

% Output the data for the experiment where we
% reduce the quasidefinite matrices from ChuLM07CAREX benchmark set

if nargin<1 || isempty(threshold)
    threshold=2;
end

% computed in analyzeChuLM07CAREX
canUseBoth = [14 15 16 19 20 28 29];
canUseFactors = [1 2 5 6 7 8 9 10 11 12 13 14 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33];

fprintf(1,'---- We use the factors B and C defining G and Q ----\n\n')

fprintf(1,'Example & $k$ & $t$ & $r$ & $\\kappa(\\mathcal{G}_{\\mathcal{I}}(X))$ & Subspace distance & $\\max|(X_{ij})|$ & $\\max|(X_{\\opt})_{ij}|$ & $\\texttt{it}$ \\\\\n\\hline\n')

for i=canUseFactors

    fprintf('%d & ',i);
    [A,G,Q,X,parout,B,R,C,Q0]=ChuLM07Carex(i);

    % these factors define X = [-Q A*; A G] = [-C*Q0C A*; A BR^{-1}B*] with
    % Q0 and R pos def

    % we use Cholesky factorizations of Q0 and R to define our factors
    % Bquad and Cquad, whose Gram matrices give G = BquadBquad* and Q =
    % Cquad*Cquad
    if i ~= 2
        cholQ0 = chol(Q0);       
    else
        cholQ0 = [3 2; 0 0];
    end
    cholR = chol(R);
    
    Bquad = B/cholR;
    Cquad = cholQ0*C;
    quad = formQuadBasis(Cquad,A,Bquad);
        % form the quadBasis representation of X;
        % some of the factors are rectangular -- this is taken care of in formQuadBasis by padding with zeros
    sym = symBasisFromQuadBasis(quad);
    [quadOpt,~,~,iterations] = optimizeQuadBasis(quad,threshold);

    U = symplecticSubspaceFromQuadBasis(quad);
    Uopt = symplecticSubspaceFromQuadBasis(quadOpt);
    symOpt = symBasisFromQuadBasis(quadOpt);

    fprintf(1, '%d & %d & %d & %1.2e & %1.2e & %5.3e & %5.3f & %d\\\\\n',...
        parout(1),parout(2),parout(3),cond(U),subspace(U,Uopt), max(max(abs(sym.X))), max(max(abs(symOpt.X))), iterations)
end

fprintf(1,'\n\n---- We use the Cholesky factors of G and Q ----\n\n')

fprintf(1,'Example & $k$ & $t$ & $r$ & $\\kappa(\\mathcal{G}_{\\mathcal{I}}(X))$ & Subspace distance & $\\max|(X_{ij})|$ & $\\max|(X_{\\opt})_{ij}|$ & \\texttt{it} \\\\\n\\hline\n')

for i=canUseBoth

    fprintf('%d & ',i);
    [A,G,Q,X,parout]=ChuLM07Carex(i);
    % these factors define X = [-Q A*; A G] with G and Q pos def

    % we use Cholesky factorizations of G and Q to define our factors
    % Bquad and Cquad
    
    Bquad = chol(G,'lower');
    Cquad = chol(Q);
    quad = formQuadBasis(Cquad,A,Bquad); % form the quadBasis representation of X
    sym = symBasisFromQuadBasis(quad);
    
    [quadOpt,~,~,iterations] = optimizeQuadBasis(quad,threshold);

    U = symplecticSubspaceFromQuadBasis(quad);
    Uopt = symplecticSubspaceFromQuadBasis(quadOpt);
    symOpt = symBasisFromQuadBasis(quadOpt);

    fprintf(1, '%d & %d & %d & %5.3e & %5.3e & %5.3f & %5.3f & %d\\\\\n',...
        parout(1),parout(2),parout(3),cond(U),subspace(U,Uopt), max(max(abs(sym.X))), max(max(abs(symOpt.X))), iterations)
end

end
