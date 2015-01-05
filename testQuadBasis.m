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

quad.X, quad.v
sym=symBasisFromQuadBasis(quad);sym.X
[quad.X, quad.v] = updateQuadBasisOut(quad.X,quad.v,2);
sym2=symBasisFromQuadBasis(quad);sym2.X
