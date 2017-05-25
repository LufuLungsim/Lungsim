%
% Automatic outlier detection for LCI data, based on relative standard
% derivation of the Data points
%
function validData  = autoDetectOutlier(x,y,validData,parameters)

        maxStdDev = parameters.Simulation.LCIbyFitMaxStdDev ;   
        s=maxStdDev+1;                                          %Be sure to make s bigger than maxStdDev to run the while loop at least once
        while s > maxStdDev && sum(validData)>5
            [params,rms,r2,handle]=multiExponentialFit(x(validData==1),y(validData==1));

            residuals=y-handle(x);                                  %Calculate how good the fit is based on the residual error
            relResiduals=abs(residuals)./y;                         %and its distribution
            s = std(relResiduals(validData==1));

            %If the data points are too scattered from the best fit, remove the worst fitting value
            if(s>maxStdDev)
               [maxDev posMaxDev]=max(relResiduals.*validData);        %Find worst fitting value (deleted values = 0)      
               validData(posMaxDev)=0;                                 %Define worst fitting value                    
            end              
        end
        
        if sum(validData)<=5
            %warndlg('autoDetectOutlier eliminated most datapoints to fit a function. Manual review of the Data is advised.');
            validData = ones(size(validData));
        end
end
