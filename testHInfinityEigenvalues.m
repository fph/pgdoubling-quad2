function testHInfinityEigenvalues

reset(RandStream.getGlobalStream);

n=6;

XH=randn(n);XH=XH+XH';
XJ=randn(n);XJ=XJ+XJ';


symh=symBasisFromSymplecticSubspace([eye(n);XH]);
symj=symBasisFromSymplecticSubspace([eye(n);XJ]);

symh=optimizeSymBasis(symh);
symj=optimizeSymBasis(symj);

UH=symplecticSubspaceFromSymBasis(symh);
UJ=symplecticSubspaceFromSymBasis(symj);

first=1:n;second=n+1:2*n;

ZH=UH(first,:);
YH=UH(second,:);

lambda=hInfinityEigenvalues(UH,UJ);

assertElementsAlmostEqual(sort(lambda),sort(eig(XH*XJ)));
