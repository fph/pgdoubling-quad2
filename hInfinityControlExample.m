function [A B1 B2 C1 C2 D11 D12 D21 D22 gammaOpt]=hInfinityControlExample(kind,parameters)
% returns several h-infinity control problems examples
%
% [A B1 B2 C1 C2 D11 D12 D21 D22 gammaOpt]=hInfinityControlExample(kind,parameters)
%
% BenBMX07Example6.1, parameters: [a]

switch kind
    case 'BenMX07Example6.1'
        if not(exist('parameters','var')) || isempty(parameters)
            parameters=1;
        end
        a=parameters;
        A=[-a 0 1 -2 1; 0 -100 0 0 0; 0 0 0 -2*a a; 0 0 0 0 1; 0 0 0 3 2];
        B1=[1;0;a;0;0];
        B2=[0; -90; 0;0;1];
        C1=[1 0 0 0 0;  0 1 0 0 0];
        C2=[0 0 1 -2 1];
        D11=zeros(2,1);
        D12=[0;1];
        D21=1;
        D22=0;
        gammaOpt=7.853923684022;
    otherwise
        error 'Wrong type';
end
