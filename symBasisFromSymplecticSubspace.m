function [symb,invcond,swaps] = symBasisFromSymplecticSubspace(U,varargin)
% symplectic basis representation of a symplectic subspace
%
% [symb,invcond,swaps]=symBasisFromSymplecticSubspace(U,options)
%
% given a (matrix whose columns span a) symplectic subspace, returns a symplectic basis representation
% i.e., symb.X and symb.v such that Pi_v * U = [I;X]
% here Pi_v = the matrix in rowSwap(?,v,'N'), and X is symmetric
%
% options:
% 'Swap': if it is set, do not attempt to optimize, but use the requested
% swap
% 'initialSwap': use the specified swap
% 'diagonalThreshold', 'offDiagonalThreshold', 'allowedInvCond'
% 'heuristic' returns heuristic and stops
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling


o = Options(varargin{:});

dt = o.get('diagonalThreshold',[]);
odt = o.get('offDiagonalThreshold',[]);

if o.isSet('heuristic')
    if o.isSet('swap') || o.isSet('initialSwap') || o.isSet('diagonalThreshold') || o.isSet('offDiagonalThreshold')
        error('PGDoubling:conflictingOptions','Conflicting options added to ''heuristics''');
    end
    [symb, invcond] = heuristicSymBasisFromSymplecticSubspace(U);
    dt = inf;
    
elseif o.isSet('swap') && ~isempty(o.look('swap'))
    if o.isSet('initialSwap') || o.isSet('diagonalThreshold') || o.isSet('offDiagonalThreshold')
        error('PGDoubling:conflictingOptions','Conflicting options added to ''swap''');
    end
    [symb,invcond] = specifiedSymBasisFromSymplecticSubspace(U,o.get('swap'));
    dt = inf; odt = inf;
elseif o.isSet('initialSwap') && ~isempty(o.get('initialSwap'))
    [symb,invcond] = specifiedSymBasisFromSymplecticSubspace(U,o.get('initialSwap'));
    if invcond<sqrt(eps(class(U)))
        %replace our guess with the heuristic
        [symb,invcond] = heuristicSymBasisFromSymplecticSubspace(U);
        if invcond < sqrt(eps(class(U)))
            warning 'PGDoubling:badSubspace' 'the provided subspace is ill-conditioned. Using O(n^3) QRPi heuristic, but it won''t help much by itself'
        else
            warning 'PGDoubling:badInitialPermutation' 'the initial guess provided is ill-conditioned, it has been replaced by the O(n^3) QRPi heuristic';
        end
    end
else %no permutation, no initial guess
    [symb,invcond] = heuristicSymBasisFromSymplecticSubspace(U);
    if invcond < sqrt(eps(class(U)))
        warning 'PGDoubling:badSubspace' 'the provided subspace is ill-conditioned. Using O(n^3) QRP heuristic, but it won''t help much by itself'
    end
end

[symb,optcond,swaps] = optimizeSymBasis(symb,dt,odt,o.get('maxSwaps',[]));
invcond = invcond * optcond;
if optcond<o.get('allowedInvCond',1e-1) %TODO: is this a good "magic value"?
    %recompute
    [symb, invcond] = specifiedSymBasisFromSymplecticSubspace(U,symb.v);
    %even if invcond is small, can't do much about it anymore
end

if invcond < sqrt(eps(class(U))) && nargout<3
    warning('PGDoubling:illConditionedMatrix', 'the matrix I am inverting has conditioning >1/sqrt(eps). This may be due to an ill-conditioned subspace');
end
