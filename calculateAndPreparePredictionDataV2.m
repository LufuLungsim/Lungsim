function [ outputModel, outputModelInterp, outputParams, LCIStopCrit] = calculateAndPreparePredictionDataV2( parameters, table, gas, predictionStructs, validData )
    

    %Extract Measurement Parameters
    nBreaths              = gas.Standard.nBreaths(1);
    startingConcentration = gas.General.cet_start(1);
    lci0025               = gas.Standard.lci(1);
    frc                   = gas.General.frcao(validData==1)*1000;
    cev                   = gas.General.cev_ds(validData==1)*1000;
    cet                   = gas.General.cet_norm(validData==1);                      
    
    %Define Headers for interpolated and absolute prediction output. 
    generalInfo = {'Breaths used', 'Breaths excluded', 'ratio of Breaths used', 'End tidal concentration'};
    header = {'LCI', 'Rel LCI Error', 'CEV', 'FRC', 'M0', 'M1', 'M2', 'R2'};
    
    offsetGeneralInfo               = length(generalInfo);
    lengthHeader                    = length(header);
    offsetHeader                    = 3;
    interpolatedPercentages         = 1.1:-0.05:0.3;
    
    %Preallocate cell size
    outputModel         = cell(size(predictionStructs, 1)+2, length(generalInfo)+lengthHeader*length(predictionStructs{1,1}));    
    outputModelInterp   = cell(length(interpolatedPercentages)+2, 1+lengthHeader*length(predictionStructs{1,1}));
    
    %Set Headers in output matrix
    outputModel(2,:)    = [generalInfo, repmat(header, 1, length(predictionStructs{1,1}))];
        
    %Set predictionStructs names to headers in output matrix
    for i = 1:length(predictionStructs{1,1})
        outputModel(1, 1+length(generalInfo)+(i-1)*lengthHeader) = {predictionStructs{1,1}(i).name};
    end
    
    %Iterate over breaths
    for i = 1:size(predictionStructs, 1)
        breath0025   = 4+i;             	%From 5 breaths to all
        cev0025      = cev(breath0025);
        %cumNetExpVol = gas.General.cumNetExpVol(breath0025);
        endTidalConcentration = cet(breath0025);
    
        %Set breath numbers and percentage
        outputModel(i+2, 1) = {breath0025};                                     %#Breaths used
        outputModel(i+2, 2) = {nBreaths-breath0025};                            %#Breaths excluded
        outputModel(i+2, 3) = {breath0025/nBreaths};                            %ratio of LCI Breaths
        outputModel(i+2, 4) = {endTidalConcentration};                          %End tidal tracer gas concentration 
        
        %Iterate over the different prediction models
        for j = 1:length(predictionStructs{1,1})
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calculate Predictions based on Model functions
            % predictionStructs(i)  1 = {'CetN2/CEV'} 
            %                       2 = {'CetN2/Num'}       
            %                       3 = {'VolExpN2/CEV'}    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                            
            predictionFunction  = predictionStructs{i,1}; 
            r2CC        = predictionFunction(j).r2(1);
            fHandleCC   = predictionFunction(j).function;
            fParamsCC   = predictionFunction(j).params;
            fCut        = @(x)fHandleCC(fParamsCC,x)-0.025*startingConcentration;	
            
            CEV2_5      = findZeroPassingPoint(fCut);
            
            if CEV2_5 < cev0025
                CEV2_5 = cev0025;
                %warning('Estimated CEV smaller than the last measured CEV, setting CEV to measured CEV %2.0f',lastMeasurementBreath );
            end
            
            %Extrapolation of FRC based on FRC calculations made for the
            %measured breaths
            coeff = polyfit(-log10(cev(1:breath0025)) , log10(frc(1:breath0025)), 1);
            FRCextrapolated = 10^coeff(2)*CEV2_5.^(-coeff(1));
            FRC = FRCextrapolated;
            %FRC         = frc(breath0025);
            
            LCI         = CEV2_5/FRC;
            LCIrelError = (LCI-lci0025)/lci0025;
            
            a = linspace(cev(breath0025), CEV2_5);
            fCet = @(x)fHandleCC(fParamsCC,x);
            cetM = [cet(1:breath0025)'                     fCet(a)];
            lciM = [cev(1:breath0025)'./frc(1:breath0025)'  a./frc(breath0025)'];
            
            [ m0, m1, m2 ] = getMomentRatios( cetM, lciM );

            saveRegion  = 1+length(generalInfo)+(j-1)*lengthHeader:1+length(generalInfo)+j*lengthHeader-1;
            outputModel(i+2, saveRegion) = num2cell([LCI, LCIrelError, CEV2_5, FRC, m0, m1, m2, r2CC]);
        end
    end
    
    %Set Annotations in output matrix
    outputModelInterp(1:2,:)                = outputModel(1:2,offsetGeneralInfo:end); 
    outputModelInterp(offsetHeader:end, 1)  = num2cell(interpolatedPercentages);
    calculatedPercentages                   = cell2mat(outputModel(offsetHeader:end, 3)); 
    
    %for each prediction Structure (Fitting method) interpolate the values calculated above 
    for i = 1:length(predictionStructs{1,1})*lengthHeader 
        valModel                    = cell2mat(outputModel(offsetHeader:end, offsetGeneralInfo+i));      
        valModel(isnan(valModel))   = 1e-16;
        valInterpModel              = interp1(calculatedPercentages, valModel, interpolatedPercentages, 'spline');
        outputModelInterp(offsetHeader:end, 1+i) = num2cell(valInterpModel);
    end
    
    %Calculate lci as defined by stop criteria
    LCIStopCrit = getLCIFromPredictionStopCriteria(outputModel, offsetGeneralInfo, lengthHeader, gas, predictionStructs, validData);
   
    %Prepare model decay time constants to save to file
    numOfParams = zeros(size(predictionFunction));
    for i = 1:length(predictionFunction)
        numOfParams(i) = length(predictionFunction(i).params);
    end
    placeHeader = cumsum([1 numOfParams(1:end)]);
    paramName = {'a', 'b', 'c', 'd', 'e'};
    outputParams = cell(length(predictionStructs)+2, sum(numOfParams));
    for i = 1:length(numOfParams)
        outputParams(1, placeHeader(i)) = {predictionStructs{1}(i).name};
        outputParams(2, placeHeader(i):placeHeader(i+1)-1) = paramName(1:numOfParams(i));
        for j = 1:length(predictionStructs)
            outputParams(j+2, placeHeader(i):placeHeader(i+1)-1) = num2cell(predictionStructs{j}(i).params);
        end                
    end
end

function val0 = findZeroPassingPoint(functionHandle)
    options = optimset('Display','off');
    for i = -5:5
        try
            val0      = fzero(functionHandle, 1*10^i, options);
            if ~isnan(val0)
                break;
            end
        catch
            
        end
    end
    
    if isnan(val0)
       val0 = 1e-16; 
    end
end