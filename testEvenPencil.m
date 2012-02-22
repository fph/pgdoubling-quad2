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

%when v=0, we should recover the Hamiltonian (up to inverting two blocks)
v=false(2*n,1);
[X,v]=evenPencil2SymBasis(AA,EE,n,m,v);
[Ah,Eh]=symBasis2HamiltonianPencil(X,v);
K=abs(jay(2*n));
assertVectorsAlmostEqual(Eh\Ah,K*H*K');

%also when v is different...
for tries=1:4 %if we try more, we stumble on a singular leading matrix
    v=logical(randi(2,2*n,1)-1);
    [X,v]=evenPencil2SymBasis(AA,EE,n,m,v);
    [Ah,Eh]=symBasis2HamiltonianPencil(X,v);
    K=abs(jay(2*n));
    assertVectorsAlmostEqual(Eh\Ah,K*H*K');
end
