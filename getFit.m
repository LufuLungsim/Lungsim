%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non-linear fit with 3 or 5 DOF, depending on the nature of the function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [params,handle,fit,rms,r2]=getFit(type,xMin,xMax,xInput,yInput)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reduction of the data to the requested interval
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x=xInput((xInput>xMin)&(xInput<xMax));
    y=yInput((xInput>xMin)&(xInput<xMax));
   
    if strcmp(type,'nonlin')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Least square fit to multi factor exponential function
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [params,rms,r2,handle]=multiExponentialFit(x,y);       
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Conventional least square fit to linear function
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [params,rms,r2]=linearFit(x,y);
        handle = @(v) params(1)+params(2)*v;
    end   
    fit=handle(xInput);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear fitting procedure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [parameters,rms,r2]=linearFit(x,y)
    
    if isempty(x) && isempty(y)
        disp('WARNING (linearFit): cannot calculate fit')
        parameters=zeros(1,2);
        rms=0;
        r2=0;
        return;
    end

    matrixA=[ones(size(x)),x];
    parameters=matrixA\y;
    fit=matrixA*parameters;
    
    if isempty(y)
        error('(linearFit): cannot calculate mean of y');
    else
        yMean=mean(y);
    end
    residuals=y-fit;
    variation=y-yMean;

    if (variation'*variation)==0
        r2=0;
    else
        r2=1-(residuals'*residuals)/(variation'*variation);
    end

    if isempty(y) || yMean==0
        rms=0;
        error('(linearFit): cannot calculate rms');
    else
        rms=sqrt((residuals'*residuals))/length(y)/abs(yMean);
    end
end