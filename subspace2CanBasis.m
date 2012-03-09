function [X,p,invcond,swaps]=subspace2CanBasis(U,varargin)
% Bounded canonical basis representation of a subspace
%
% [X,p,invcond]=subspace2CanBasis(U,options);
%
% input: U subspace, p initial permutation guess (may be empty, in this
% case the procedure will use a O(n^3) heuristic based on QRP)
%
% output:
% X,p: such that U(p(1:end),:) and [I;X] span the same subspace
%
% Options:
% 'permutation': if this parameter exists, used its value "p" as the
% permutation and skips optimization overall (i.e., works exactly as
% subspace2SpecifiedCanBasis)
% 'initialPermutation', initial permutation guess  (may be empty, uses QRP
%  heuristic if so or if the guess is bad)
% 'threshold', threshold to use
% 'allowedInvCond', recomputes the basis if it notices an inverse condition
% number (as in Matlab's linsolve) below a certain tolerance

o=matgic.Options(varargin{:});

if o.isSet('permutation')
    [X,p,invcond]=subspace2SpecifiedCanBasis(U,o.get('permutation'));
else
    if o.isSet('initialPermutation') && ~isempty(o.get('initialPermutation'))
        [X,p,invcond]=subspace2SpecifiedCanBasis(U,o.get('initialPermutation'));
        if invcond<sqrt(eps(class(U)))
            %replace our guess with the heuristic
            [X,p,invcond]=subspace2SpecifiedCanBasis(U,canBasisHeuristic(U));
            if invcond<sqrt(eps(class(U)))
                warning 'cbrpack:badSubspace' 'subspace2CanBasis: the provided subspace is ill-conditioned. Using O(n^3) QRP heuristic, but it won''t help much by itself'
            else
                warning 'cbrpack:badInitialPermutation' 'subspace2CanBasis: the initial guess provided is ill-conditioned, it has been replaced by the O(n^3) QRP heuristic';
            end
        end
    else %no permutation, no initial guess
        [X,p,invcond]=subspace2SpecifiedCanBasis(U,canBasisHeuristic(U));
        if invcond<sqrt(eps(class(U)))
            warning 'cbrpack:badSubspace' 'subspace2CanBasis: the provided subspace is ill-conditioned. Using O(n^3) QRP heuristic, but it won''t help much by itself'
        end
    end
    [X,p,optcond,swaps]=optimizeCanBasis(X,p,o.get('threshold',[]),o.get('maxSwaps',[]));
    if optcond<o.get('allowedInvCond',1e-1) %TODO: is this a good "magic value"?
        %recompute
        [X,p,invcond]=subspace2SpecifiedCanBasis(U,p);
        %even if invcond is small, can't do much about it anymore
    end
end



