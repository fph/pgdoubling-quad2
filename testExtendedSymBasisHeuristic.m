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
        
    EE=[U(1:n,:)' zeros(n+m,m)];
    AA=[(jay(n)*U(n+1:2*n,:))' U(2*n+1:2*n+m,:)'];
%    size(AA),size(EE),size(v)
    
    S=warning('off','cbrpack:notSymplectic');
    [X,v,invcond]=evenPencil2SymBasis(AA,EE,n/2,m,vh);
    [Xp,vp,invcondp]=evenPencil2SymBasis(AA,EE,n/2,m,vhp);
    warning(S);
    
    assertElementsAlmostEqual(invcond,1);
    assertElementsAlmostEqual(invcondp,1);

    assertEqual(logical(vhp(:)),logical(v(:)));
    assertEqual(logical(vh(:)),logical(v(:)));
    
end
