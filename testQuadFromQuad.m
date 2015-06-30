function testQuadFromQuad()

% test updateQuadBasisOut vs. quadBasisFromQuadBasis in this special case
quad.X = [1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16];
quad.v = [true true true false];
sym = symBasisFromQuadBasis(quad);

newv = [true false true false];

[firstQuad.X, firstQuad.v] = updateQuadBasisOut(quad.X,quad.v,2);

firstSym = symBasisFromQuadBasis(firstQuad);
symForward = symBasisFromSymBasis(sym, firstSym.v);

assertEqual(firstQuad.v, newv);
assertElementsAlmostEqual(symForward.X, firstSym.X);

secondQuad = quadBasisFromQuadBasis(quad,newv);

secondSym = symBasisFromQuadBasis(secondQuad);
symForward = symBasisFromSymBasis(sym, secondSym.v);

assertEqual(secondQuad.v, newv);
assertElementsAlmostEqual(symForward.X, secondSym.X);

% randomized tests updateQuadBasisOut vs. quadBasisFromQuadBasis in this special case
reset(RandStream.getGlobalStream);
n = 5;
for i = 1:50
    quad.X = rand(n) + 1i*rand(n);
    quad.v = logical(randi([0,1],n,1));
    out = randi(n);
    quad.v(out) = true;
    newv = quad.v;
    newv(out) = false;
    sym = symBasisFromQuadBasis(quad);
    
    [firstQuad.X, firstQuad.v] = updateQuadBasisOut(quad.X,quad.v,out);
    
    firstSym = symBasisFromQuadBasis(firstQuad);
    symForward = symBasisFromSymBasis(sym, firstSym.v);
    
    assertEqual(firstQuad.v, newv);
    assertElementsAlmostEqual(symForward.X, firstSym.X);

    secondQuad = quadBasisFromQuadBasis(quad,newv);
    
    secondSym = symBasisFromQuadBasis(secondQuad);
    symForward = symBasisFromSymBasis(sym, secondSym.v);
    
    assertEqual(secondQuad.v, newv);
    assertElementsAlmostEqual(symForward.X, secondSym.X);
end

% test updateQuadBasisIn vs. quadBasisFromQuadBasis in this special case
quad.X = [1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16];
quad.v = [true false false false];
sym = symBasisFromQuadBasis(quad);

newv = [true false true false];

[firstQuad.X, firstQuad.v] = updateQuadBasisIn(quad.X,quad.v,3);

firstSym = symBasisFromQuadBasis(firstQuad);
symForward = symBasisFromSymBasis(sym, firstSym.v);

assertEqual(firstQuad.v, newv);
assertElementsAlmostEqual(symForward.X, firstSym.X);

secondQuad = quadBasisFromQuadBasis(quad,newv);

secondSym = symBasisFromQuadBasis(secondQuad);
symForward = symBasisFromSymBasis(sym, secondSym.v);

assertEqual(secondQuad.v, newv);
assertElementsAlmostEqual(symForward.X, secondSym.X);

% randomized tests updateQuadBasisIn vs. quadBasisFromQuadBasis in this special case
reset(RandStream.getGlobalStream);
n = 5;
for i = 1:50
    quad.X = rand(n) + 1i*rand(n);
    quad.v = logical(randi([0,1],n,1));
    in = randi(n);
    quad.v(in) = false;
    newv = quad.v;
    newv(in) = true;
    sym = symBasisFromQuadBasis(quad);
    
    [firstQuad.X, firstQuad.v] = updateQuadBasisIn(quad.X,quad.v,in);
    
    firstSym = symBasisFromQuadBasis(firstQuad);
    symForward = symBasisFromSymBasis(sym, firstSym.v);
    
    assertEqual(firstQuad.v, newv);
    assertElementsAlmostEqual(symForward.X, firstSym.X);

    secondQuad = quadBasisFromQuadBasis(quad,newv);
    
    secondSym = symBasisFromQuadBasis(secondQuad);
    symForward = symBasisFromSymBasis(sym, secondSym.v);
    
    assertEqual(secondQuad.v, newv);
    assertElementsAlmostEqual(symForward.X, secondSym.X);
end

% test updateQuadBasisInOut vs. quadBasisFromQuadBasis in this special case
quad.X = [1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16]+1i*rand(4);
quad.v = [true true false false];
sym = symBasisFromQuadBasis(quad);

newv = [true false true false];

[firstQuad.X, firstQuad.v] = updateQuadBasisInOut(quad.X,quad.v,3,2);

firstSym = symBasisFromQuadBasis(firstQuad);
symForward = symBasisFromSymBasis(sym, firstSym.v);

assertEqual(firstQuad.v, newv);
assertElementsAlmostEqual(symForward.X, firstSym.X);

secondQuad = quadBasisFromQuadBasis(quad,newv);

secondSym = symBasisFromQuadBasis(secondQuad);
symForward = symBasisFromSymBasis(sym, secondSym.v);

assertEqual(secondQuad.v, newv);
assertElementsAlmostEqual(symForward.X, secondSym.X);

% randomized test updateQuadBasisInOut vs. quadBasisFromQuadBasis in this special case
reset(RandStream.getGlobalStream);
n = 5;
for i = 1:50
    quad.X = rand(n)+ 1i*rand(n);
    quad.v = logical(randi([0,1],n,1));
    in = randi(n);
    out = randi(n);
    if in == out
        continue;
    end
    quad.v(in) = false;
    quad.v(out) = true;
    newv = quad.v;
    newv(in) = true;
    newv(out) = false;    
    sym = symBasisFromQuadBasis(quad);
    
    [firstQuad.X, firstQuad.v] = updateQuadBasisInOut(quad.X,quad.v,in,out);
    
    firstSym = symBasisFromQuadBasis(firstQuad);
    symForward = symBasisFromSymBasis(sym, firstSym.v);
    
    assertEqual(firstQuad.v, newv);
    assertElementsAlmostEqual(symForward.X, firstSym.X);

    secondQuad = quadBasisFromQuadBasis(quad,newv);
    
    secondSym = symBasisFromQuadBasis(secondQuad);
    symForward = symBasisFromSymBasis(sym, secondSym.v);
    
    assertEqual(secondQuad.v, newv);
    assertElementsAlmostEqual(symForward.X, secondSym.X);
end

% randomized test for quadBasisFromQuadBasis Case 1 with |K| > 1.

reset(RandStream.getGlobalStream);
n = 10;
for i = 1:50
    quad.X = rand(n) + 1i*rand(n);
    quad.v = logical(randi([0,1],n,1));
    newv = logical( quad.v.*logical(randi([0,1],n,1)) ); %\mathcal{J} \subseteq  \mathcal{I}
    
    sym = symBasisFromQuadBasis(quad);
      
    quad = quadBasisFromQuadBasis(quad,newv);
    sym2 = symBasisFromQuadBasis(quad);
    
    symForward = symBasisFromSymBasis(sym, sym2.v);
    
    assertEqual(quad.v, newv);
    assertElementsAlmostEqual(symForward.X, sym2.X);
end

% randomized test for quadBasisFromQuadBasis Case 2 with |K| > 1.

reset(RandStream.getGlobalStream);
n = 10;
for i = 1:50
    quad.X = rand(n) + 1i*rand(n);
    newv = logical(randi([0,1],n,1));
    quad.v = logical( quad.v.*logical(randi([0,1],n,1)) ); %\mathcal{J} \supseteq  \mathcal{I}
    
    sym = symBasisFromQuadBasis(quad);
      
    quad = quadBasisFromQuadBasis(quad,newv);
    sym2 = symBasisFromQuadBasis(quad);
    
    symForward = symBasisFromSymBasis(sym, sym2.v);
    
    assertEqual(quad.v, newv);
    assertElementsAlmostEqual(symForward.X, sym2.X);
end

% randomized test for quadBasisFromQuadBasis Case 3 with |K| > 2.

reset(RandStream.getGlobalStream);
n = 10; k = 0;
for i = 1:50
    quad.X = rand(n) + 1i*rand(n);    
    
    while true
        quad.v = logical(randi([0,1],n,1));
        w = logical( quad.v.*logical(randi([0,1],n,1)) );
        if any(w ~= false(n,1)) && any(quad.v ~= true(n,1)), break, end
    end
    
    newv = ~w;
      
    sym = symBasisFromQuadBasis(quad);
      
    quad = quadBasisFromQuadBasis(quad,newv);
    sym2 = symBasisFromQuadBasis(quad);
    
    symForward = symBasisFromSymBasis(sym, sym2.v);
    
    assertEqual(quad.v, newv);
    assertElementsAlmostEqual(symForward.X, sym2.X);
end

