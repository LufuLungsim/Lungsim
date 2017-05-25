%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculation of a pure fractional delay according to the Thiran transfer 
% function, see e.g., http://www.mathworks.com/help/toolbox/control/ref/thiran.html
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 1.0, 15. April 2011
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [b]=filterLagrange(delay,N)
%
% implements a lagrange interpolation filter of order N for fractional delay d
%
    D=max(delay);               % maximal scaled delay
    Nfull=N*floor(D/N);         % separation of full orders of digital delays
    Neff=Nfull+N;               % effective filter order
    
    len=length(delay);
    
    b=zeros(Neff+1,len);        % nominator coeffients
    
    for i=1:len
        D=delay(i);
        Nfull=N*floor(D/N);         % separation of full orders of digital delays
        Deff=D-Nfull;               % fractional delay
        for k=0:N                   % calculation of lagrange interpolation coefficients
            nTerms=Deff-[0:N];
            nTerms(k+1)=1;
            dTerms=k-[0:N];
            dTerms(k+1)=1;
            b(Nfull+k+1,i)=prod(nTerms./dTerms);
        end
    end
end
