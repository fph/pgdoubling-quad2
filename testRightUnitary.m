function testRightUnitary

randn('seed',0);

m = 5;
n = 6;

for tries=1:100

  A = rand(m,n);
  rows = logical(randi([0 1],m,1));
  k = sum(rows);
  cols = [ones(1,n-k) zeros(1,k)];
  cols = cols(randperm(n));

  [H,Anew] = rightUnitary(A,rows,cols);

  assertVectorsAlmostEqual(H*H',eye(size(H)));
  assertVectorsAlmostEqual(A*H,Anew);
  assertVectorsAlmostEqual(Anew(rows==1,cols==1),zeros(k,n-k));

end