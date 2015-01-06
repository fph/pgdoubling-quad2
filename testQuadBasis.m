function testQuadBasis()

% test symBasisFromQuadBasis
quad.X = [1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16];
quad.v = [true true true false];
sym = symBasisFromQuadBasis(quad);
assertEqual(quad.v, sym.v);
Y =[-107  -122  -137    13
    -122  -140  -158    14
    -137  -158  -179    15
      13    14    15   256];
assertElementsAlmostEqual(sym.X,Y);

% test updateQuadBasisOut
sym = symBasisFromQuadBasis(quad);
[quad.X, quad.v] = updateQuadBasisOut(quad.X,quad.v,2);
assertEqual(quad.v, [true false true false]);
sym2 = symBasisFromQuadBasis(quad);
symForward = symBasisFromSymBasis(sym, sym2.v);
assertElementsAlmostEqual(symForward.X, sym2.X);

%test updateQuadBasisIn
quad.X = [1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16];
quad.v = [true false false false];
sym = symBasisFromQuadBasis(quad);
[quad.X, quad.v] = updateQuadBasisIn(quad.X,quad.v,3);
assertEqual(quad.v, [true false true false]);
sym2 = symBasisFromQuadBasis(quad);
symForward = symBasisFromSymBasis(sym, sym2.v);
assertElementsAlmostEqual(symForward.X, sym2.X);

% randomized tests updateQuadBasisOut
reset(RandStream.getGlobalStream);
n = 5;
for i = 1:50
    quad.X = rand(n) + 1i*rand(n);
    quad.v = logical(randi([0,1],n,1));
    out = randi(n);
    quad.v(out) = true;
    sym = symBasisFromQuadBasis(quad);
    [quad.X, quad.v] = updateQuadBasisOut(quad.X,quad.v,out);
    sym2 = symBasisFromQuadBasis(quad);
    symForward = symBasisFromSymBasis(sym, sym2.v);
    assertElementsAlmostEqual(symForward.X, sym2.X);
end

% randomized tests updateQuadBasiIn
reset(RandStream.getGlobalStream);
n = 5;
for i = 1:50
    quad.X = rand(n) + 1i*rand(n);
    quad.v = logical(randi([0,1],n,1));
    in = randi(n);
    quad.v(in) = false;
    sym = symBasisFromQuadBasis(quad);
    [quad.X, quad.v] = updateQuadBasisIn(quad.X,quad.v,in);
    sym2 = symBasisFromQuadBasis(quad);
    symForward = symBasisFromSymBasis(sym, sym2.v);
    assertElementsAlmostEqual(symForward.X, sym2.X);
end

