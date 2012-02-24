function testLeftDual

n=10;

X=randn(n);
p=randperm(2*n);

U=canBasis2Subspace(X,p);

[Etilde,Atilde]=leftDual(X,p);

assertVectorsAlmostEqual([Etilde Atilde]*jay(2*n)*U,zeros(n));
