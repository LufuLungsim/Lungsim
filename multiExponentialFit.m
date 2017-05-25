%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multi-exponential fit based on the paper "Fast and Accurate Fitting and Filtering of Noisy
% Exponentials in Legendre Space"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The non-linear fit is a double exponential ansatz function that is found by fitting
% first to a single exponential with constant term, used as initialization.
% Stability is satisfacorily.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

function [p1,rms,r2,handle]=multiExponentialFit(x,y)

    if isempty(x) && isempty(y)
        fprintf('WARNING (multiExponentialFit): cannot calculate fit')
        p1=[];
        rms=0;
        r2=0;
        handle=@(v) 0;
        return;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fitting strategy: trying fits for 3 exponential components first. If
    % lambda structure is not phyological, trying with p=2 and then with
    % p=1 exponentials
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    expFit = ExpFit;
    m=3+floor(length(y)/5);
    p=3;
    [lambda,coefficients,zFit ] = expFit.do(x, y, m, p);
    if lambda(1)>0 || lambda(2)>0 || lambda(3)>0 || norm(imag(lambda))>0 || coefficients(2)<0 || coefficients(3)<0 || coefficients(4)<0
        p=2;
        [lambda,coefficients,zFit ] = expFit.do(x, y, m, p);
        if lambda(1) > 0 || lambda(2) > 0 || norm(imag(lambda)) > 0 || coefficients(2)<0 || coefficients(3)<0
            fprintf('WARNING (multiExponentialFit): only single exponential fit possible\n')
            p=1;
            [lambda,coefficients,zFit ] = expFit.do(x, y, m, p);
        end
    end

    handle=@(v) exp(v*[0;lambda]')*coefficients;
    p1=[coefficients(1),reshape([coefficients(2:end),lambda]',1,2*p)]';

    %     figure(99)
    %     plot(x,y,x*2,handle(x*2))
    %
    yMean=mean(y);
    residuals=y-handle(x);
    variation=y-yMean;
           
    if (variation'*variation)==0
        r2=0;
    else
        r2=1-(residuals'*residuals)/(variation'*variation);
    end
    
    if isempty(y) || yMean==0
        rms=0;
        error('(multiExponentialFit): cannot calculate relative rms');
    else
        rms=sqrt((residuals'*residuals))/length(y)/yMean;
    end
end