function testExtendedSymBasisHeuristic

reset(RandStream.getGlobalStream);

n=6; m=2;

for i=1:50
    
    %creates a matrix with only one "good" permutation
    U=[[eye(n) zeros(n,m)];zeros(n,n+m);[zeros(m,n) eye(m)]]+1e-10*randn(2*n+m,n+m);
    v=logical(randi(2,n,1)-1);
    U(1:2*n,:)=rowSwap(U(1:2*n,:),v,'N');
    
    [vh invcondprev]=extendedSymBasisHeuristicOld(U);
    [vhp invcondprevp]=extendedSymBasisHeuristicPaper(U);
    
    first = 1:n/2;
    second = n/2+1:n;
    third = n+1:3*n/2;
    fourth = 3*n/2+1:2*n;
    
    EE=[-U(second,:)' -U(third,:)' zeros(n+m,m)];
    AA=[U(first,:)' U(fourth,:)' U(2*n+1:2*n+m,:)'];
%    size(AA),size(EE),size(v)
    
    [sym,invcond]=symBasisFromEvenPencil(AA,EE,n/2,m,vh);
    [symp,invcondp]=symBasisFromEvenPencil(AA,EE,n/2,m,vhp);
    
    assertElementsAlmostEqual(invcond,1);
    assertElementsAlmostEqual(invcondp,1);

    assertEqual(logical(vhp(:)),logical(symp.v(:)));
    assertEqual(logical(vh(:)),logical(sym.v(:)));
    
end
