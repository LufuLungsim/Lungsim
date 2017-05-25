%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculation of a butterworth filter coefficients 
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 2.0, 25. August 2011
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a,b]=filterButterworth(fCrit,N,tSample)
%
% implements a discret butterworth filter of order N and bandwidth fCrit
%
    wCrit=fCrit*2*pi;           % critical angular frequency
    
    d=zeros(N+1,1);             % denominator coeffients
    n=1;                        % nominator coeffients
    
    if N==1                     % order 1
        d=[1/wCrit,1];
    elseif N==2                 % order 2
        d=[1/wCrit^2,sqrt(2)/wCrit,1];
    elseif N==3                 % order 3
        d=[1/wCrit^3,2/wCrit^2,2/wCrit,1];
    else
        fprintf('Butterworth filter of order %d not yet implemented!\n',N);
    end
    
    [a,b]=filterBilinear(n,d,0,tSample);
end
    
