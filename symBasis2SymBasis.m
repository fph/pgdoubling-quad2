function [X,targetv]=symBasis2SymBasis(X,v,targetv)
% updates a symBasis to a new specified "target" v
%
% [X,targetv]=symBasis2SymBasis(X,v,targetv)
%
% equivalent to symplecticSubspace2SymBasis(symBasis2SymplecticSubspace(X,v),'swap',targetv);
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

%crappy implementation for now...
[X, targetv]=symplecticSubspace2SymBasis(symBasis2SymplecticSubspace(X,v),'swap',targetv);
