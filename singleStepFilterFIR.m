%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Implementation of FIR filter using optional initialization (steady state)
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 1.0, 19. Nov. 2014
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y]=singleStepFilterFIR(b,x,j)
%
% implements a FIR filter for output element j, initialized to the steady state assuming x(k)=x(1) for k<1.
%
lengthB=length(b);
% starting from x(j), the elements until x(j-lengthB+1) are used in xi. 
% If vector x is not long enough, xi is padded with its first entry x(1)
xi=[x(j:-1:1+max(0,j-lengthB));ones(max(0,lengthB-j),1)*x(1)];
y=xi'*b(:);

