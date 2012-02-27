function lambda=hInfinityEigenvalues(UH,UJ)
% computes eigenvalues of X_H*X_J for a h-Infinity control problem
%
% given bases of the subspaces UH=im [X_H;I], UJ=im [X_J;I], computes the
% eigenvalues of (X_HX_J)

[twon n]=size(UH);
first=1:n;second=n+1:2*n;

ZH=UH(first,:);
YH=UH(second,:);

ZJ=UJ(first,:);
YJ=UJ(second,:);

U=[YH';ZJ'];
[X,p]=subspace2CanBasis(U);
[Ytilde Ztilde]=leftDual(X,p);

Ytilde=Ytilde';
Ztilde=Ztilde';

%assertVectorsAlmostEqual(YH*Ztilde,ZJ*Ytilde)

lambda=eig(ZH*Ztilde,YJ*Ytilde);
