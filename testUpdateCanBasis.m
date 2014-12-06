function testUpdateCanBasis

reset(RandStream.getGlobalStream);

U=randn(10,6);
U(2,1)=10^12; %so that we test a "difficult" case

ip=[3 4 1 9 7 10 6 5 2 8];

can=canBasisFromSubspace(U,'permutation',ip);
assertEqual(ip,can.p);

[newcan.X newcan.p]=updateCanBasis(can.X,can.p,[3 4,1],[6,5,1]);
newcan.origin='subspace';
assertEqual(newcan.p,[6 4 1 9 8 2 3 5 10 7]);

U2=subspaceFromCanBasis(newcan);
assertElementsAlmostEqual(subspace(U2,U),0);

U=[111 12; 113 14; 115 116; 17 118; 119 110]; %for this matrix, the QRP heuristic leads to a max(max(abs(X)))>1
can=canBasisFromSubspace(U,'threshold',10);
x=max(max(abs(can.X)));
assert(x<10);
assert(x>1.1);

can=canBasisFromSubspace(U,'threshold',1.2);
x=max(max(abs(can.X)));
assert(x<1.2);
assert(x>1.1);

can=canBasisFromSubspace(U,'threshold',1.1);
x=max(max(abs(can.X)));
assert(x<1.1);

for k=1:100
    U=randn(40,25);
    can=canBasisFromSubspace(U,'threshold',1.2);
    U2=subspaceFromCanBasis(can);
    assertElementsAlmostEqual(subspace(U,U2),0);
end
