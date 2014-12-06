function testDoublingStep

reset(RandStream.getGlobalStream);

n=10;

for tries=1:10
    sym.X=randn(n);sym.X=(sym.X+sym.X')/2;
    sym.v=logical(randi(2,n,1)-1);
    sym.origin='symplecticPencil';
    [L,U]=symplecticPencilFromSymBasis(sym);
    sym2=doublingStep(sym);
    M=U\L;
    sym3=symBasisFromSymplecticPencil(M^2,eye(n));
    
    [LL,UU]=symplecticPencilFromSymBasis(sym2);
    [LLL,UUU]=symplecticPencilFromSymBasis(sym3);
    
    assertElementsAlmostEqual(subspace([LL,UU]',[LLL,UUU]'),0);
end
