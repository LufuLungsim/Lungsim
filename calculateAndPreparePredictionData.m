function [ outputModel,  outputMeasurement, outputModelInterp, outputMeasurementInterp, LCIStopCrit] = calculateAndPreparePredictionData( parameters, table, gas, predictionStructs, validData )
    

    %Extract Measurement Parameters
    nBreaths              = gas.Standard.nBreaths(1);
    startingConcentration = gas.General.cet_start(1);
    lci0025               = gas.Standard.lci(1);
    cev                   = gas.General.cev_ds(validData==1);
    cet                   = gas.General.cet_norm(validData==1);                      
    
    %Define Headers
    generalInfo = {'Breaths used', 'Breaths excluded', 'ratio of Breaths used', 'End tidal concentration'};
    header = {'LCI', 'Rel LCI Error', 'Breath 2_5', 'CEV', 'VN2Exp', 'FRC', 'R2 %/CEV', 'R2 %/#', 'R2 Vol/CEV', 'Estimated Error'};
    
    offsetGeneralInfo               = length(generalInfo);
    lengthHeader                    = length(header);
    offsetHeader                    = 3;
    interpolatedPercentages         = 1.1:-0.05:0.3;
    
    %Preallocate cell size
    outputModel         = cell(size(predictionStructs, 1)+2, length(generalInfo)+lengthHeader*length(predictionStructs{1,1}));
    outputMeasurement   = cell(size(predictionStructs, 1)+2, length(generalInfo)+lengthHeader*length(predictionStructs{1,1}));
    
    outputModelInterp               = cell(length(interpolatedPercentages)+2, 1+lengthHeader*length(predictionStructs{1,1}));
    outputMeasurementInterp         = cell(length(interpolatedPercentages)+2, 1+lengthHeader*length(predictionStructs{1,1}));
    
    %Set Headers in output matrix
    outputModel(2,:)        = [generalInfo, repmat(header, 1, length(predictionStructs{1,1}))];
    outputMeasurement(2,:)  = [generalInfo, repmat(header, 1, length(predictionStructs{1,1}))];
        
    %Set predictionStructs names to headers in output matrix
    for i = 1:length(predictionStructs{1,1})
        outputModel(1, 1+length(generalInfo)+(i-1)*lengthHeader) = {predictionStructs{1,1}(i).name};
        outputMeasurement(1, 1+length(generalInfo)+(i-1)*lengthHeader) = {predictionStructs{1,1}(i).name};
    end

    for i = 1:size(predictionStructs, 1)
        breath0025   = 4+i;             	%From 5 breaths to all
        cev0025      = cev(breath0025);
        cumNetExpVol = gas.General.cumNetExpVol(breath0025);
        endTidalConcentration = cet(breath0025);
    
        %Set breath numbers and percentage
        outputModel(i+2, 1) = {breath0025};                                     %#Breaths used
        outputModel(i+2, 2) = {nBreaths-breath0025};                            %#Breaths excluded
        outputModel(i+2, 3) = {breath0025/nBreaths};                            %ratio of LCI Breaths
        outputModel(i+2, 4) = {endTidalConcentration};                          %End tidal tracer gas concentration 
        outputMeasurement(i+2, 1) = outputModel(i+2, 1);                        %#Breaths used
        outputMeasurement(i+2, 2) = outputModel(i+2, 2);                        %#Breaths excluded
        outputMeasurement(i+2, 3) = outputModel(i+2, 3);                        %ratio of Breaths
        outputMeasurement(i+2, 4) = outputModel(i+2, 4);                        %End tidal tracer gas concentration         
        
        for j = 1:length(predictionStructs{1,1})
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calculate Predictions based on Model functions
            % predictionStructs(i,  1 = {'CetN2/CEV'} 
            %                       2 = {'CetN2/Num'}       
            %                       3 = {'VolExpN2/CEV'}    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                            
            predictionFunction  = predictionStructs{i,1}; 
            r2CC        = predictionFunction(j).r2(1);
            fHandleCC   = predictionFunction(j).function;
            fParamsCC   = predictionFunction(j).params;
            fCut        = @(x)fHandleCC(fParamsCC,x)-0.025*startingConcentration;
            CEV2_5      = findZeroPassingPoint(fCut);
                                                                  
            predictionFunction  = predictionStructs{i,2}; 
            r2CN        = predictionFunction(j).r2(1);
            fHandleCN   = predictionFunction(j).function;
            fParamsCN   = predictionFunction(j).params;
            fCut        = @(x)fHandleCN(fParamsCN,x)-0.025*startingConcentration;
            Breath2_5   = findZeroPassingPoint(fCut);
            
            predictionFunction  = predictionStructs{i,3}; 
            r2VN        = predictionFunction(j).r2(1);
            fHandleVN   = predictionFunction(j).function;
            fParamsVN   = predictionFunction(j).params;
            
            Breaths     = linspace(0, Breath2_5);
            VN2Exp      = trapz(Breaths, fHandleVN(fParamsVN, Breaths));
          
            FRC         = VN2Exp/(startingConcentration-0.025*startingConcentration);

            LCI         = CEV2_5/FRC;
            LCIrelError = (LCI-lci0025)/lci0025;
            
            saveRegion  = 1+length(generalInfo)+(j-1)*lengthHeader:1+length(generalInfo)+(j)*lengthHeader-1;
            outputModel(i+2, saveRegion) = num2cell([LCI, LCIrelError, Breath2_5, CEV2_5, VN2Exp, FRC, r2CC, r2CN, r2VN, 3-r2CC-r2CN-r2VN]);
            
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calculate Predictions based on Model functions and
            % Measuremet data
            % predictionStructs(i,  1 = {'CetN2/CEV'} 
            %                       2 = {'CetN2/Num'} 
            %                       3 = {'VolExpN22/Num'}    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            
            predictionFunction           = predictionStructs{i,1}; 
            r2CC        = predictionFunction(j).r2(1);
            fHandleCC   = predictionFunction(j).function;
            fParamsCC   = predictionFunction(j).params;
            fCut        = @(x)fHandleCC(fParamsCC,x)-0.025*startingConcentration;
            CEV2_5      = findZeroPassingPoint(fCut);
            
            %Make sure the estimated CEV is larger than the last measured
            %one
            if CEV2_5 < cev0025
                CEV2_5 = cev0025;
                %warning('Estimated CEV smaller than the last measured CEV, setting CEV to measured CEV %2.0f',lastMeasurementBreath );
            end
            
            predictionFunction           = predictionStructs{i,2}; 
            r2CN        = predictionFunction(j).r2(1);
            fHandleCN   = predictionFunction(j).function;
            fParamsCN   = predictionFunction(j).params;
            fCut        = @(x)fHandleCN(fParamsCN,x)-0.025*startingConcentration;
            Breath2_5   = findZeroPassingPoint(fCut);
            
            if Breath2_5 < breath0025
                Breath2_5 = breath0025;
                %warning('Last breath larger than estimated number of Breaths needed to reach N2 2.5 %2.0f',lastMeasurementBreath );
            end
            
            predictionFunction           = predictionStructs{i,3}; 
            r2VN        = predictionFunction(j).r2(1);
            fHandleVN   = predictionFunction(j).function;
            fParamsVN   = predictionFunction(j).params;
            
            Breaths     = linspace(breath0025, Breath2_5);
            VN2ExpModel = trapz(Breaths, fHandleVN(fParamsVN, Breaths));
            VN2ExpMeas  = cumNetExpVol*1000;
            VN2Exp      = VN2ExpModel+VN2ExpMeas;
            
            FRC         = VN2Exp/(startingConcentration-0.025*startingConcentration);

            LCI         = CEV2_5/FRC;
            LCIrelError = (LCI-lci0025)/lci0025;
            
            saveRegion  = 1+length(generalInfo)+(j-1)*lengthHeader:1+length(generalInfo)+j*lengthHeader-1;
            outputMeasurement(i+2, saveRegion) = num2cell([LCI, LCIrelError, Breath2_5, CEV2_5, VN2Exp, FRC, r2CC, r2CN, r2VN, 3-r2CC-r2CN-r2VN]);
        end
    end
    
    %Set Annotations in output matrix
    outputModelInterp(1:2,:)        = outputModel(1:2,offsetGeneralInfo:end); 
    outputMeasurementInterp(1:2,:)  = outputMeasurement(1:2,offsetGeneralInfo:end);
    
    outputModelInterp(offsetHeader:end, 1) = num2cell(interpolatedPercentages);
    outputMeasurementInterp(offsetHeader:end, 1) = num2cell(interpolatedPercentages);
    
    calculatedPercentages           = cell2mat(outputModel(offsetHeader:end, 3)); 
    
    for i = 1:length(predictionStructs{1,1})*lengthHeader         %for each prediction Structure (Fitting method) interpolate the values calculated above 
        valModel                    = cell2mat(outputModel(offsetHeader:end, offsetGeneralInfo+i));      
        valModel(isnan(valModel))   = 1e-16;
        valInterpModel              = interp1(calculatedPercentages, valModel, interpolatedPercentages, 'spline');
        outputModelInterp(offsetHeader:end, 1+i) = num2cell(valInterpModel);

        valMeas                 = cell2mat(outputMeasurement(offsetHeader:end, offsetGeneralInfo+i));
        valMeas(isnan(valMeas)) = 1e-16;
        valInterpMeasurement    = interp1(calculatedPercentages, valMeas, interpolatedPercentages, 'spline');
        outputMeasurementInterp(offsetHeader:end, 1+i) = num2cell(valInterpMeasurement);
    end
    
    LCIStopCrit = getLCIFromPredictionStopCriteria(outputMeasurement, offsetGeneralInfo, lengthHeader, gas, predictionStructs);
   
end

function val0 = findZeroPassingPoint(functionHandle)
    options = optimset('Display','off');
    for i = -3:3
        val0      = fzero(functionHandle, 1*10^i, options);
        if ~isnan(val0)
            break;
        end
    end
    
    if isnan(val0)
       val0 = 1e-16; 
    end
end