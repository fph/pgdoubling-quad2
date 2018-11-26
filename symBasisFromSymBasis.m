function symb=symBasisFromSymBasis(symb,targetv)
% updates a symBasis to a new specified "target" v
%
% sym = symBasis2SymBasis(symb,targetv)
%
% equivalent to symplecticSubspace2SymBasis(symBasis2SymplecticSubspace(symb),'swap',targetv);
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

range = 1:length(targetv);
inout = range(symb.v~=targetv);
[symb.X, symb.v] = updateSymBasis(symb.X, symb.v, inout);
