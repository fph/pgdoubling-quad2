function [X,targetv]=symBasis2symBasis(X,v,targetv)
% updates a symBasis to a new specified "target" v
%
% [X,targetv]=symBasis2symBasis(X,v,targetv)
%
% equivalent to symplecticSubspace2SymBasis(symBasis2SymplecticSubspace(X,v),'swap',targetv);

%crappy implementation for now...
[X, targetv]=symplecticSubspace2SymBasis(symBasis2SymplecticSubspace(X,v),'swap',targetv);
