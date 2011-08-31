function [X,p,swaps]=subspace2CanBasis(U,varargin)
% Canonical basis representation of a subspace
%
% [X,p,swaps]=subspace2CanBasis(U);
%
% output: 
% X,p: such that U(p(1:end),:) and [I;X] span the same subspace
% swaps: number of "swaps" performed during the algorithm
% (each swap is essentially a rank-1 update to X and costs O(n^2))
%
% optional arguments:
%
% 'initialPermutation',p: a guess for the starting permutation p
% If omitted, it is replaced by a O(n^3) guess based on QRP
%
% 'threshold',T (default 2): the maximum allowed magnitude of the elements in X. 
% May be infinity; in this case, the starting permutation (or the heuristic) is used
%
% 'recompute',k (default 2): recomputes X and v at the end from the initial data, for added stability, if k or more swaps are needed
% May be 0 (always recompute) or inf (never)
%
% 'maxSwaps',M (default 10*length(U)): never uses more than M swaps, warns if threshold criterion is not met afterwards
%

o=matgic.Options(varargin{:});

maxSwaps=o.get('maxSwaps',10*length(U));
threshold=o.get('threshold',2);
if threshold<1+sqrt(eps(U))
    error('cbrpack:thresholdTooSmall','you can only hope to enforce thresholds T=1+sqrt(eps) or larger');
end

k=o.get('recompute',2);
if(o.isSet('initialPermutation'))
    p=o.get('initialPermutation');
else
    p=canBasisHeuristic(U);
end

[m n]=size(U);
assert(length(p)==m);

%matlab crappy boilerplate to catch warnings...
saved=warning('query','MATLAB:singularMatrix');
warning off MATLAB:singularMatrix;
lastwarn('');

X=U(p(n+1:m),:)/U(p(1:n),:);

%...and to rethrow them
[lastmsg,lastid]=lastwarn;
if(strcmp(lastid,'MATLAB:singularMatrix'))
    warning('cbrpack:badInitialPermutation', 'bad initial permutation provided - computing a new one using the heuristic. All fine, but this costs O(n^3)');
    p=canBasisHeuristic(U);
    X=U(p(n+1:m),:)/U(p(1:n),:);
end
warning(saved);

%optimization loop

swaps=0;
while(swaps<maxSwaps)
    [maxvec, maxis]=max(abs(X));
    [maxval maxj]=max(maxvec);
    maxi=maxis(maxj);
    % the three lines above compute maxi,maxj=argmax(abs(X(i,j)))

    if maxval<threshold
        break;
    end

    [X,p]=updateCanBasis(X,p,maxi,maxj);
    swaps=swaps+1;
end
if swaps==maxSwaps
    warning('cbrpack:stagnated','failed to produce a X with elements below the required threshold (obtained:%d, required:%d). Try running with a larger threshold.',maxval,threshold);
end

if swaps>=k
    X=U(p(n+1:m),:)/U(p(1:n),:); %recompute X based on p
end
