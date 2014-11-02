function testSymBasisFromSymplecticSubspace
reset(RandStream.getGlobalStream);

U=[1:6]';
assertExceptionThrown(@() rowSwap(U,[1 0 0 0],'N'),'PGDoubling:wrongPermutationLength');
assertExceptionThrown(@() rowSwap(U,[1 0 0],'S'),'PGDoubling:wrongTransposedness');

assertEqual(rowSwap(U,[0 0 0],'N'),[1 2 3 4 5 6]');
assertEqual(rowSwap(U,[0 0 0],'T'),[1 2 3 4 5 6]');

assertEqual(rowSwap(U,[1 1 1],'N'),[4 5 6 -1 -2 -3]');

assertEqual(rowSwap(U,[1 0 0],'N'),[4 2 3 -1 5 6]');
assertEqual(rowSwap(U,[1 0 0],'T'),[-4 2 3 1 5 6]');

U=[2;1];
symb=symBasisFromSymplecticSubspace(U);
assertEqual(symb.origin,'symplecticSubspace');

assertExceptionThrown(@() symBasisFromSymplecticSubspace([1]),'PGDoubling:oddSize');

U=randn(10,4);
assertExceptionThrown(@() symBasisFromSymplecticSubspace(U),'PGDoubling:oddSize');

n=8;
for k=1:10
    U=randomLagrangianSubspace(2*n);
    v1=randi([0,1],n,1);
    v2=randi([0,1],n,1);
    symb1=symBasisFromSymplecticSubspace(U,'swap',v1);
    assertEqual(symb1.v,v1);
    symb2=symBasisFromSymplecticSubspace(U,'swap',v2);
    assertEqual(symb2.v,v2);
    range=1:n; inout=range(v1~=v2);
    [symb12.X,symb12.v]=updateSymBasis(symb1.X,symb1.v,inout);
    assertEqual(symb12.v,symb2.v);
    assertElementsAlmostEqual(symb12.X,symb2.X);
end

n=8;
for k=1:10
    U=randomLagrangianSubspace(2*n);
    symb=symBasisFromSymplecticSubspace(U);
    U2=symplecticSubspaceFromSymBasis(symb);
    assertElementsAlmostEqual(subspace(U2,U),0);
end

n=8;
for S=[100 5 4 3 2 1.1 1.01 1.0001]
    U=randomLagrangianSubspace(2*n);
    symb=symBasisFromSymplecticSubspace(U,'diagonalThreshold',S);
    T=sqrt(5+S^2);
    assert(all(all(abs(symb.X)<=T)));
    U2=symplecticSubspaceFromSymBasis(symb);
    assertElementsAlmostEqual(subspace(U2,U),0);
end
