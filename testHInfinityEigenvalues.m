function testHInfinityEigenvalues

reset(RandStream.getGlobalStream);

n=6;

XH=randn(n);XH=XH+XH';
XJ=randn(n);XJ=XJ+XJ';


[XHP,vh]=symplecticSubspace2SymBasis([eye(n);XH]);
[XJP,vj]=symplecticSubspace2SymBasis([eye(n);XJ]);

[XHP,vh]=optimizeSymBasis(XHP,vh);
[XJP,vj]=optimizeSymBasis(XJP,vj);

UH=symBasis2SymplecticSubspace(XHP,vh);
UJ=symBasis2SymplecticSubspace(XJP,vj);

first=1:n;second=n+1:2*n;

ZH=UH(first,:);
YH=UH(second,:);

lambda=hInfinityEigenvalues(UH,UJ);

assertElementsAlmostEqual(sort(lambda),sort(eig(XH*XJ)));
