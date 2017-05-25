%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Implementation of FIR filter with non-constant coefficients
% uses optional initialization (steady state)
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 1.3, 05. Mai 2011
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y]=filterFIR(b,x)
%
% implements a FIR filter for optionally non-constant coefficients, initialized to the steady state 
% assuming x(k)=x(1) for k<1.
%
    sizeB=size(b);
    
    y=zeros(size(x));
    for i=1:length(x)
        if sizeB(2)==1	% checks if constant filter coefficients present
            k=1;
        else
            k=i;
        end
        if i<sizeB(1)
            xi=[x(i:-1:1);ones(sizeB(1)-i,1)*x(1)];
        else
            xi=x(i:-1:i-sizeB(1)+1);
        end
        y(i)=xi'*b(:,k);
    end
end
