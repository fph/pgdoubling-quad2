function [can,invcond,swaps]=canBasisFromSubspace(U,varargin)
% Bounded canonical basis representation of a subspace
%
% [can,invcond,swaps]=subspace2CanBasis(U,options);
%
% input: U subspace
%
% output:
% can: such that U(can.p(1:end),:) and [I;can.X] span the same subspace
% invcond: inverse condition number
% swaps: number of swaps performed after the initial permutation (specified
%        or heuristic)
%
% Options:
% 'permutation': if this parameter exists, use its value "p" as the
% permutation and skips optimization overall
% 'initialPermutation', initial permutation guess  (may be empty, uses O(n^3)
%  QRP heuristic if so or if the guess is bad)
% 'threshold', threshold to use
% 'allowedInvCond', recomputes the basis if it notices an inverse condition
% number (as in Matlab's linsolve) below a certain tolerance
% 'heuristic': returns the QRP heuristic without further optimization
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

o = Options(varargin{:});

threshold = o.get('threshold',[]);

if o.isSet('heuristic')
    if o.isSet('permutation') || o.isSet('initialPermutation') || o.isSet('threshold')
        error('PGDoubling:conflictingOptions','Conflicting options added to ''heuristics''');
    end
    [can, invcond] = heuristicCanBasisFromSubspace(U);
    threshold = inf;
elseif o.isSet('permutation')
    if o.isSet('initialPermutation') || o.isSet('threshold')
        error('PGDoubling:conflictingOptions','Conflicting options added to ''permutation''');
    end
    [can,invcond]=specifiedCanBasisFromSubspace(U,o.get('permutation'));
    threshold = inf;
elseif o.isSet('initialPermutation')
    [can, invcond] = specifiedCanBasisFromSubspace(U,o.get('initialPermutation'));
    if invcond < sqrt(eps(class(U)))
        %replace our guess with the heuristic
        [can, invcond] = heuristicCanBasisFromSubspace(U);
        if invcond < sqrt(eps(class(U)))
            warning 'PGDoubling:badSubspace' 'canBasisFromSubspace: the provided subspace is ill-conditioned. Using O(n^3) QRP heuristic, but it won''t help much by itself'
        else
            warning 'PGDoubling:badInitialPermutation' 'canBasisFromSubspace: the initial guess provided is ill-conditioned, it has been replaced by the O(n^3) QRP heuristic';
        end
    end
else
    [can, invcond] = heuristicCanBasisFromSubspace(U);
    if invcond < sqrt(eps(class(U)))
        warning 'PGDoubling:badSubspace' 'canBasisFromSubspace: the provided subspace is ill-conditioned. Using O(n^3) QRP heuristic, but it won''t help much by itself'
    end
end

[can, optcond, swaps]=optimizeCanBasis(can, threshold, o.get('maxSwaps',[]));
if optcond < o.get('allowedInvCond',1e-1) %TODO: is this a good "magic value"?
    %recompute
    [can, invcond] = specifiedCanBasisFromSubspace(U, p);
    %even if invcond is small, can't do much about it anymore
end

