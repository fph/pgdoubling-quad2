function [X,v,invcond,swaps]=symplecticSubspace2SymBasis(U,varargin)
% symplectic basis representation of a symplectic subspace
%
% [X,v,invcond,swaps]=symplecticSubspace2SymBasis(U,options)
%
% given a (matrix whose columns span a) symplectic subspace, returns a symplectic basis representation
% i.e., X and v such that Pi_v * U = [I;X]
% here Pi_v = the matrix in rowSwap(?,v,'N'), and X is symmetric
%
% options:
% 'Swap': if it is set, do not attempt to optimize, but use the requested
% swap
% 'initialSwap': use the specified swap
% 'diagonalThreshold', 'offDiagonalThreshold', 'allowedInvCond'
% 


o=matgic.Options(varargin{:});

if o.isSet('swap')
    [X,v,invcond]=symplecticSubspace2SpecifiedSymBasis(U,o.get('swap'));
    swaps=0;
else
    if o.isSet('initialSwap') && ~isempty(o.get('initialSwap'))
        [X,v,invcond]=symplecticSubspace2SpecifiedSymBasis(U,o.get('initialSwap'));
        if invcond<sqrt(eps(class(U)))
            %replace our guess with the heuristic
            [X,v,invcond]=symplecticSubspace2HeuristicSymBasis(U);
            if invcond<sqrt(eps(class(U)))
                warning 'cbrpack:badSubspace' 'subspace2SymBasis: the provided subspace is ill-conditioned. Using O(n^3) QRP heuristic, but it won''t help much by itself'
            else
                warning 'cbrpack:badInitialPermutation' 'subspace2SymBasis: the initial guess provided is ill-conditioned, it has been replaced by the O(n^3) QRP heuristic';
            end
        end
    else %no permutation, no initial guess
        [X,v,invcond]=symplecticSubspace2HeuristicSymBasis(U);
        if invcond<sqrt(eps(class(U)))
            warning 'cbrpack:badSubspace' 'subspace2SymBasis: the provided subspace is ill-conditioned. Using O(n^3) QRP heuristic, but it won''t help much by itself'
        end
    end
    [X,v,optcond,swaps]=optimizeSymBasis(X,v,o.get('diagonalThreshold',[]),o.get('offDiagonalThreshold',[]),o.get('maxSwaps',[]));
    if optcond<o.get('allowedInvCond',1e-1) %TODO: is this a good "magic value"?
        %recompute
        [X,v,invcond]=symplecticSubspace2SpecifiedSymBasis(U,v);
        %even if invcond is small, can't do much about it anymore
    end
end

if invcond<sqrt(eps(class(U))) && nargout<3
    warning('cbrpack:illConditionedMatrix', 'symplecticSubspace2SymBasis: the matrix I am inverting has conditioning >1/sqrt(eps). This may be due to an ill-conditioned subspace or to a bad initial guess --- consider using the initial value heuristic instead');
end
