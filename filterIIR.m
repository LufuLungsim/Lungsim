%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Implementation of MATLAB filter with initialization (steady state)
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 1.0, 12. April 2011
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y,zf]=filterIIR(b,a,x,zi)
%
% implements a filter that is initialized to the the steady state 
% assuming x(k)=x(1) for k<1.
%

    b=b/a(1);                       % scale coefficient vectors to make a(1)=1
    a=a/a(1);
    
    len=max(length(a),length(b));   % determine maximal degree of the polynomials
    
    a(length(a)+1:len,1)=0;         % pad the vectors with zeros
    b(length(b)+1:len,1)=0;
    
    if isempty(zi)                % initialisation, if zi==[]
        rhs=[b,-a]*[1;sum(b)/sum(a)];
        S=diag(ones(1,len-1))+diag(ones(1,len-2)*-1,1);
        if length(rhs)<2
            zi=[];
        else
            zi=(S\rhs(2:end))*x(1);
        end
    end
    
    [y,zf]=filter(b,a,x,zi);        % filtering with MATLAB function
end