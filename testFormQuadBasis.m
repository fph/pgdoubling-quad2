function testFormQuadBasis

randn('seed',0);

m = 6;
n = 5;

% test that for square B and C formQuadBasis just saves them in X
C = rand(n,n);
A = rand(m,n);
B = rand(m,m);

X = [C zeros(n,m); A B];
quadF = formQuadBasis(C,A,B);
  
assertEqual(quadF.v,[true(1,n) false(1,m)]);
assertElementsAlmostEqual(quadF.X,X);

% test for rectangular factors

for i = 1:50
    C = rand(n,2*n);
    A = rand(2*m,2*n);
    B = rand(2*m,m);
    quadF = formQuadBasis(C,A,B);
  
    assertEqual(quadF.v,[true(1,2*n) false(1,2*m)]);
    
    X = quadF.X;
    v = quadF.v;
    Csquare = X(v,v);
    Bsquare = X(~v,~v);
    
    assertEqual(Csquare,[zeros(n,2*n);C])
    assertEqual(Bsquare,[B zeros(2*m,m)])
    assertElementsAlmostEqual(Csquare'*Csquare,C'*C)
    assertElementsAlmostEqual(Bsquare*Bsquare',B*B')
    assertEqual(X(~v,v),A)
end
