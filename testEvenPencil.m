function testEvenPencil

reset(RandStream.getGlobalStream);

[A,G,Q,X,parout,B,R]=carex(3);

[n m]=size(B);

H=hamiltonian(A,G,Q);
[AA,EE]=evenPencil(A,B,Q,R);

%since S=0, the systems should be the same
v1=eig(AA,EE);
v2=eig(H);
assertElementsAlmostEqual(sort(abs(v1(isfinite(v1)))),sort(abs(v2)));

%We should recover the Hamiltonian (up to inverting two blocks)
for tries=1:2 %if we try more, we stumble on a singular leading matrix
    v=logical(randi(2,2*n,1)-1);
    sym=symBasisFromEvenPencil(AA,EE,n,m,v);
    [Ah,Eh]=hamiltonianPencilFromSymBasis(sym);
    K=abs(jay(2*n));
    assertVectorsAlmostEqual(Eh\Ah,K*H*K');
end
