%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Appliying bilinear transformation to continous transfer function 
% Basis of this approach see: http://www.mathworks.com/help/toolbox/signal/bilinear.html
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 1.1, 27. April 2011
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a,b]=filterBilinear(n,d,wrapFraction,Tsample)
%
% implements first order dynamic system by bilinear transformation
%
    fs=1/Tsample;                   % calculation of prewarping frequency
    if wrapFraction ~= 0
        fp=2*pi*fs/2*wrapFraction;	% scaling of maximal frequency
        fs=fp/tan(fp/fs/2);
    else
        fs=2*fs;
    end
    
    p=roots(d);                     % roots of denominator
    z=roots(n);                     % roots of nominator
    k=n(1)/d(1);                    % gain
    
    pd = (1+p/fs)./(1-p/fs);        % bilinear transformation
    zd = (1+z/fs)./(1-z/fs);
    kd = real(k*prod(fs-z)./prod(fs-p));
    
    z0 = ones(1,length(p))*(-1);	% insertion of zero at -1 (to ensure correct working in z-space)
    z0(find(zd)) = zd;
    zd = z0;
    
    b=poly(zd)'*kd;                 % nominiator (as column vector)
    a=poly(pd)';                    % dnominiator (as column vector)
end
