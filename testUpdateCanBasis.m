function testUpdateCanBasis

reset(RandStream.getDefaultStream);

U=randn(10,6);
U(2,1)=10^12; %so that we test a "difficult" case

ip=[3 4 1 9 7 10 6 5 2 8];

[X,p]=subspace2CanBasis(U,'threshold',inf,'initialPermutation',ip);
assertEqual(ip,p);

[Xn,pn]=updateCanBasis(X,p,[3 4,1],[6,5,1]);
assertEqual(pn,[6 4 1 9 8 2 3 5 10 7]);

U2=canBasis2Subspace(Xn,pn);
assertElementsAlmostEqual(subspace(U2,U),0);
