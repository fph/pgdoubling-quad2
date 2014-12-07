function testInverseFreeSignStep
reset(RandStream.getGlobalStream)
for trie = 1:100
    n=6;
    sym.X=randn(n);sym.X=sym.X+sym.X';sym.v=logical(randi(2,n,1)-1);
    sym.origin='hamiltonianPencil';
    
    [A,E]=hamiltonianPencilFromSymBasis(sym);
    H=E\A;
    Hnew=1/2*(H+inv(H));
    
    [symnew,w,swaps1,swaps2,res]=inverseFreeSignStep(sym);
    
    [Anew,Enew]=hamiltonianPencilFromSymBasis(symnew);
    assertElementsAlmostEqual(0,subspace([Anew';Enew'],[Hnew';eye(size(Hnew))]));
end