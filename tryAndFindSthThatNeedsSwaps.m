1;
n=200;
for i=1:100000
    S=randomLagrangianSubspace(n);
    [X,v,invcond,swaps]=symplecticSubspace2SymBasis(S,'diagonalThreshold',1+1e-7,'offDiagonalThreshold',sqrt(2)+1e-7);
    if swaps>0
        save 'asd' S;
        'found!'
        break;
    end
end
