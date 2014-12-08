n=2;
reset(RandStream.getGlobalStream)
G = randn(n);G=G*G';
H = randn(n);H=H*H';
A = randn(n);
Ham = hamiltonian(A,G,H);

SH = Ham*jay(2*n);
sym.X = SH;
sym.v=[0 0 1 1];
symp = symBasisFromSymBasis(sym,[0,0,0,0]);
symp.X

sym = symBasisFromHamiltonianPencil(Ham,eye(2*n),'swap',[0 0 0 0]);
