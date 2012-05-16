function [A B1 B2 C1 C2 D11 D12 D21 D22 gammaOpt]=hInfinityControlExample(kind,parameters)
% returns several h-infinity control problems examples
%
% [A B1 B2 C1 C2 D11 D12 D21 D22 gammaOpt]=hInfinityControlExample(kind,parameters)
%
% BenBMX07Example6.1, parameters: [a]
% DrcexcExample2.4, parameters: [m c k]
% BenMX07Example3.1: parameters [alpha beta delta epsilon1 epsilon2 eta]
% BenMX07Example6.2 or 6.3: as the previous one, with parameters hardcoded to the
%   values in the paper
%
% (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

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
    case 'DrcexcExample2.4'
        % Uncertainty model of the Mass/Damper/Spring system
        
        if not(exist('parameters','var')) || isempty(parameters)
            parameters=[3 1 2];
        end
        m=parameters(1);
        c=parameters(2);
        k=parameters(3);
        
        pm = 0.4;
        pc = 0.2;
        pk = 0.3;
        %
        A = [
            0 1
            -k/m -c/m];
        B1 = [ 0 0 0
            -pm -pc/m -pk/m];
        B2 = [ 0
            1/m];
        C1 = [-k/m -c/m
            0 c
            k 0 ];
        C2 = [ 1 0 ];
        D11 = [ -pm -pc/m -pk/m
            0 0 0
            0 0 0];
        D12=[1/m
            0
            0];
        D21 = [0 0 0];
        D22 = 0;
        gammaOpt=nan;
    case 'BenMX07Example6.4'
        if not(exist('parameters','var')) || isempty(parameters)
            parameters=3;
        end
        alpha=parameters(1);

        A=blkdiag(2,-1);
        B1=[0 1; 0 1];
        B2=[-1; -2];
        C1=eye(2);
        D11=blkdiag(alpha,-1);
        D12=[0;1];
        C2=[4 -2];
        D21=[0 1];
        D22=0;
        gammaOpt=alpha;
    case 'BenMX07Example3.1'
        if not(exist('parameters','var')) || isempty(parameters)
            error 'You should specify the parameter values'
        end
        alpha=parameters(1);
        beta=parameters(2);
        delta=parameters(3);
        epsilon1=parameters(4);
        epsilon2=parameters(5);
        eta=parameters(6);
        
        A=-eye(2);
        B1=blkdiag(epsilon1,epsilon2);
        B2=[1;1];
        C1=blkdiag(alpha,beta);
        C2=[delta eta];
        D11=blkdiag(1/2,1/2);
        D12=[0;1];
        D21=D12';
        D22=0;
        gammaOpt=nan;
    case 'BenMX07Example6.2'
        [A B1 B2 C1 C2 D11 D12 D21 D22]=hInfinityControlExample('BenMX07Example3.1',[1 1 1 0 0 1]);
        gammaOpt=1/2;
    case 'BenMX07Example6.3'
        [A B1 B2 C1 C2 D11 D12 D21 D22]=hInfinityControlExample('BenMX07Example3.1',[1 1 1 0 1 1]);
        gammaOpt=0.8062257748299;
    otherwise
        error 'Wrong type';
end
