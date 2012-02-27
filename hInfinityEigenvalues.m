function lambda=hInfinityEigenvalues(UH,UJ)
% computes eigenvalues of X_H*X_J for a h-Infinity control problem
%
% given bases of the subspaces UH=im [X_H;I], UJ=im [X_J;I], computes the
% eigenvalues of (X_HX_J)

[twon n]=size(UH);
first=1:n;second=n+1:2*n;

%these are as follows and not the opposite following (25) in BenBMX07. I
%suspect they tweaked the definition of the even pencils in order to
%maintain the same order.
YH=UH(first,:);
ZH=UH(second,:);

YJ=UJ(first,:);
ZJ=UJ(second,:);

U=[YH';ZJ'];
[X,p]=subspace2CanBasis(U);
[Ytilde Ztilde]=leftDual(X,p);

Ytilde=Ytilde';
Ztilde=Ztilde';

%assertVectorsAlmostEqual(YH*Ztilde,ZJ*Ytilde)

ZH*Ztilde, YJ*Ytilde

lambda=eig(ZH*Ztilde,YJ*Ytilde);
