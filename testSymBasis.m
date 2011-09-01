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

U=randn(10,5).*exp(randn(10,5));
[X,v]=symplecticSubspace2SymBasis(U);
U2=symBasis2SymplecticSubspace(X,v);
assertElementsAlmostEqual(subspace(U2,U),0);


