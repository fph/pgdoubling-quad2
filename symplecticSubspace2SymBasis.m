function [X,v,swaps]=symplecticSubspace2SymBasis(U,varargin)
% symplectic basis representation of a symplectic subspace
%
% [X,v,swaps]=symplecticSubspace2SymBasis(U,...options...)
%
% given a (matrix whose columns span a) symplectic subspace, returns a symplectic basis representation
% i.e., X and v such that Pi_v * U = [I;X]
% here Pi_v = the matrix in rowSwap(?,v,'N'), and X is symmetric
%
% swaps is the number of swaps performed during the optimization algorithm.
% Swapping two rows at the same time counts as 2.
%
% optional arguments:
%
% 'initialRowSwap',p: a guess for the starting row swap v
% If omitted, it is replaced by a O(n^3) guess based on QRP
%
% 'diagonalThreshold',S: the maximum allowed magnitude of diag(X). 
%
% 'offDiagonalThreshold', T: the maximum allowed magnitude of the other elements of X.
% The theory requires S>=1, T>=sqrt(1+S^2); you'd better be far from those
% values to be ok in critical cases. Defaults are (S=2, T=sqrt(5+S^2)).
% Note that 'diagonalThreshold',inf is allowed (don't optimize at all, just return the starting basis).
%
% 'recompute',k (default 2): recomputes X and v at the end from the initial data, for added stability, if k or more swaps are needed
% May be 0 (always recompute) or inf (never)
%
% 'maxSwaps',M (default 10*length(U)): never uses more than M swaps, warns if threshold criterion is not met afterwards

o=matgic.Options(varargin{:});

maxSwaps=o.get('maxSwaps',10*length(U));
diagonalThreshold=o.get('diagonalThreshold',2);
offDiagonalThreshold=o.get('offDiagonalThreshold',sqrt(5+diagonalThreshold^2));

if diagonalThreshold<1+sqrt(eps(class(U)))
    error('cbrpack:thresholdTooSmall','you can only hope to enforce thresholds S=1+sqrt(eps) or larger');
end

if offDiagonalThreshold<sqrt(1+diagonalThreshold^2)+sqrt(eps(class(U)))
    error('cbrpack:thresholdTooSmall','you can only hope to enforce thresholds T=sqrt(1+S^2) or larger');
end

k=o.get('recompute',2);
if(o.isSet('initialRowSwap'))
    v=o.get('initialRowSwap');
else
    v=symBasisHeuristic(U);
end

[m n]=size(U);
if m~=2*n
    error('cbrpack:oddSize','symplectic matrices should be 2nxn');
end

S=rowSwap(U,v,'N');

%matlab crappy boilerplate to catch warnings...
saved=warning('query','MATLAB:singularMatrix');
warning off MATLAB:singularMatrix;
lastwarn('');

X=S(n+1:end,:)/S(1:n,:);

%...and to rethrow them
[lastmsg,lastid]=lastwarn;
if(strcmp(lastid,'MATLAB:singularMatrix'))
    warning('cbrpack:badInitialPermutation', 'bad initial permutation provided - computing a new one using the heuristic. All fine, but this costs O(n^3)');
    v=symBasisHeuristic(U);
    S=rowSwap(U,v,'N');
    X=S(n+1:end,:)/S(1:n,:);
end
warning(saved);

%optimization
swaps=0;
while(swaps<maxSwaps)
    [maxvec, maxis]=max(abs(X));
    [maxval maxj]=max(maxvec);
    maxi=maxis(maxj);
    % the three lines above compute maxi,maxj=argmax(abs(X(i,j)))

    [maxdiag maxdiagPos]=max(diag(X));
    if maxdiag>diagonalThreshold
        [X,v]=updateSymBasis(X,v,maxdiagPos);
        swaps=swaps+1;
    elseif maxval>offDiagonalThreshold
        [X,v]=updateSymBasis(X,v,[maxi maxj]);
        swaps=swaps+2;
    else
        break;
    end    
end

if swaps==maxSwaps
    warning('cbrpack:stagnated','failed to produce a X with elements below the required threshold (obtained:%d, required:%d). Try running with a larger threshold.',maxval,threshold);
end

if swaps>=k
    U=rowSwap(U,v,'N');
    X=U(n+1:end,:)/S(1:n,:);
end

%final check
if norm(X-X','fro')/norm(X) > sqrt(eps)
    warning('cbrpack:notSymplectic','the resulting matrix is numerically very far from Hermitian --- was your starting subspace symplectic?');
end

%symmetrize
X=(X+X')/2;
