function testSymBasis
reset(RandStream.getDefaultStream);

U=[1:6]';
assertExceptionThrown(@() rowSwap(U,[1 0 0 0],'N'),'cbrpack:wrongPermutationLength');
assertExceptionThrown(@() rowSwap(U,[1 0 0],'S'),'cbrpack:wrongTransposedness');

assertEqual(rowSwap(U,[0 0 0],'N'),[1 2 3 4 5 6]');
assertEqual(rowSwap(U,[0 0 0],'T'),[1 2 3 4 5 6]');

assertEqual(rowSwap(U,[1 1 1],'N'),[4 5 6 -1 -2 -3]');

assertEqual(rowSwap(U,[1 0 0],'N'),[4 2 3 -1 5 6]');
assertEqual(rowSwap(U,[1 0 0],'T'),[-4 2 3 1 5 6]');

U=[2;1];
[X,v]=symplecticSubspace2SymBasis(U);

assertExceptionThrown(@() symplecticSubspace2SymBasis([1]),'cbrpack:oddSize');

U=randn(10,4);
assertExceptionThrown(@() symplecticSubspace2SymBasis(U),'cbrpack:oddSize');

n=8;
for k=1:10
    U=randomLagrangianSubspace(2*n);
    v1=randi([0,1],n,1);
    v2=randi([0,1],n,1);
    [X1,v]=symplecticSubspace2SymBasis(U,'diagonalThreshold',inf,'initialRowSwap',v1);
    assertEqual(v,v1);
    [X2,v]=symplecticSubspace2SymBasis(U,'diagonalThreshold',inf,'initialRowSwap',v2);
    assertEqual(v,v2);
    range=[1:n];inout=range(v1~=v2);
    [X,v]=updateSymBasis(X1,v1,inout);
    assertEqual(v,v2);
    assertElementsAlmostEqual(X,X2);
end

U=randomLagrangianSubspace(2*n);
[X,v]=symplecticSubspace2SymBasis(U);
U2=symBasis2SymplecticSubspace(X,v);
assertElementsAlmostEqual(subspace(U2,U),0);
