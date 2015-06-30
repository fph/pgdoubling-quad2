function testLeftUnitary

randn('seed',0);

m = 6;
n = 5;

for tries=1:100

  A = rand(m,n);
  cols = logical(randi([0 1],1,n));
  k = sum(cols);
  rows = [ones(1,m-k) zeros(1,k)];
  rows = rows(randperm(m));

  [H,Anew] = leftUnitary(A,rows,cols);

  assertVectorsAlmostEqual(H*H',eye(size(H)));
  assertVectorsAlmostEqual(H*A,Anew);
  assertVectorsAlmostEqual(Anew(rows==1,cols==1),zeros(m-k,k));

end