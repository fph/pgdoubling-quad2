function [A,G,Q,X, parout, B, R, C, Q0]=ChuLM07Carex(n)
%generates the nth experiment in ChuLM07, n=1..33

experiments={
    1,[],
    2,[],
    3,[],
    4,[],
    5,[],
    6,[],
    7,[1],
    7,[1e-6],
    8,[1],
    8,[1e-8],
    9,[1],
    9,[1e6],
    9,[1e-6],
    10,[1],
    10,[1e-5],
    10,[1e-7],
    11,[1],
    11,[0],
    12,[1],
    12,[1e6],
    13,[1],
    13,[1e-6],
    14,[1],
    14,[1e-6],
    15,[39],
    15,[119],
    15,[199],
    16,[8],
    16,[64],
    17,[21,1,1],
    17,[21,100,100],
    18,[], %n=100 is the default
    19,[]
};

assertEqual(size(experiments),[33,2]);
assert(n>=1 && n<=33);

if(isempty(experiments{n,2}))
    [A,G,Q,X parout B R C Q0]=carex(experiments{n,1}); %need to switch due to a bug in carex(19,[])
else
    [A,G,Q,X parout B R C Q0]=carex(experiments{n,1},experiments{n,2});
end

