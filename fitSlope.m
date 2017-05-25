%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Fitting procedure in two variants: 
%       lin:        least square to line
%       normlin:    least square to line of normed values 
%       nonlin:     fitting to sum of exponentials (with adjustment of decay parameter)
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 1.0, 14. Mai 2014
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fitParameters,fit,rms,rSquared,errorValue]=fitSlope(type,xMin,xMax,x,y)

    errorValue = 0;
    
    indices=find(1==(x>xMin).*(x<xMax));
    if length(indices)<4
        fitParameters=[];
        fit=[];
        rms=0;
        rSquared=0;
        errorValue = 1;
        return;
    end
    yData=y(indices);
    
    if strcmp(type,'lin') || strcmp(type,'normlin')
        matrixA=[ones(size(x)),x];
        matrixDataA=matrixA(indices,:);
        fitParameters=matrixDataA\yData;
        fit=matrixA*fitParameters;
    else
        error('(fitSlope): unknown type=%s found\n', type);
    end
    
    if isempty(yData)
        error('(fitSlope): cannot calculate mean of yData');
    else
        yDataMean=mean(yData);
    end
    residuals=yData-fit(indices);
    variation=yData-yDataMean;

    if (variation'*variation)==0
        rSquared=0;
    else
        rSquared=1-(residuals'*residuals)/(variation'*variation);
    end

    if isempty(yData) || yDataMean==0
        rms=0;
        error('(fitSlope): cannot calculate rms');
    else
        rms=sqrt((residuals'*residuals))/length(yData)/abs(yDataMean);
    end
end

