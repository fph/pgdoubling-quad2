function testLeftDual

n=10;

can.X=randn(n);
can.p=randperm(2*n);
can.origin='subspace';

U=subspaceFromCanBasis(can);

[Etilde,Atilde]=leftDual(can);

assertVectorsAlmostEqual([Etilde Atilde]*jay(2*n)*U,zeros(n));
