function testInverseQuasidef

% 2x2 case; compare against the know formula 

C = 2;
A = 1;
B = 1;
[N,K,L,~] = inverseQuasidef(C,A,B);
T = [-N*N' K; K' L'*L]; % factors nonunique, inverse unique
trueInv = 1/(abs(B*C)^2 + abs(A)^2)*[-B^2 A; A C^2];
assertElementsAlmostEqual(T, trueInv);

C = 0;
B = 0;
A = 1;
[N,K,L,~] = inverseQuasidef(C,A,B);
T = [-N*N' K; K' L'*L];
trueInv = [0 1; 1 0];
assertElementsAlmostEqual(T, trueInv);

% randomized tests
reset(RandStream.getGlobalStream);
m = 3;
n = 5;
for i = 1 : 50
   C = rand(n,n);
   B = rand(m,m);
   A = rand(m,n);
   P = [-C'*C A'; A B*B'];
   Pinv = inv(P);
   [N,K,L,~] = inverseQuasidef(C,A,B);
   T = [-N*N' K; K' L'*L];
   assertElementsAlmostEqual(T, Pinv);
end

% randomized tests; C and B low rank
reset(RandStream.getGlobalStream);
m = 4;
n = 5;
for i = 1 : 50
   C = rand(n,2)*rand(2,n);
   B = rand(m,3)*rand(3,m);
   A = rand(m,n);
   P = [-C'*C A'; A B*B'];
   Pinv = inv(P);
   [N,K,L,~] = inverseQuasidef(C,A,B);
   T = [-N*N' K; K' L'*L];
   assertElementsAlmostEqual(T, Pinv);
end