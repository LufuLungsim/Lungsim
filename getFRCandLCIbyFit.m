% The LCI by Fit routine was implemented to have a more robust method to 
% determine the LCI. The advantage of the by Fit method is that the LCI is 
% not defined by a single breath but by the combination of breaths. 
% As defined in the consensus Paper  the LCI is defined as when the 
% concentration of the tracer gas is less than 1/40 of the original 
% concentration for three consecutive expirations. From the three consecutive 
% breaths the first is taken to define the LCI. With the by Fit method a 
% linear interpolation between the last breath above the 1/40 threshold and 
% the first below is taken, and the intersection with the exact concentration 
% point is found by linear interpolation to find the CEV and FRC

function gas = getFRCandLCIbyFit(validData, table, gas, compatibility, parameters, Washin)
    
    if Washin 
        mode = 'Washin';
    else
        mode = 'Washout';
    end

    disp(' ');
    disp(['LCI calculation by fitting Cet(norm) and FRC for ' gas.name ' for ' mode]);

    momentRatio = parameters.Simulation.MomentRatio;
    linFit      = parameters.Simulation.LinearByFit;
    expFit      = parameters.Simulation.ExponentialByFit;
    
    %Predefine the size of some vectors
    gas.ByFit.cet_norm = zeros(size(table.criticalEndRatios));
    gas.ByFit.lci      = zeros(size(table.criticalEndRatios));
    gas.ByFit.frc      = zeros(size(table.criticalEndRatios));
    gas.ByFit.duration = zeros(size(table.criticalEndRatios));
    gas.ByFit.nBreaths = zeros(size(table.criticalEndRatios));

    validCet_norm     = gas.General.cet_norm(validData==1);
    validCEV_DS       = gas.General.cev_ds(validData==1); 
    
    %Get the exponential functions of the measurements 
    if expFit
        handleCet       = getExpFit(validCEV_DS, validCet_norm); 
        handleFRC       = getExpFit(validCEV_DS, gas.General.frcao(validData==1));
        handleDuration  = getExpFit(validCEV_DS, cumsum(table.breathDuration(validData==1)));
    end
    
    %Calculate FRC & LCI, searching for the last breath above the
    %treshold and the first below, interpolate linearly between those two points and set FCR and LCI at
    %the intersection
    for i = 1:length(table.criticalEndRatios) 
                
        CetCrit = table.criticalEndRatios(i);
        
        %find the last point above the critial end ratio while the three
        %following points must be below            
        indices           = find(validCet_norm<CetCrit);            
        pointBelowCetCrit = getLCIIndex(indices,compatibility);        
        
        if pointBelowCetCrit > 1
            pointAboveCetCrit=pointBelowCetCrit-1;
            CEVmax = validCEV_DS(pointAboveCetCrit+1);
      
            %Get the linar functions of the points surrounding the by Fit point 
            if linFit 
                handleCet       = getLinFit(pointAboveCetCrit, validCEV_DS, validCet_norm); 
                handleFRC       = getLinFit(pointAboveCetCrit, validCEV_DS, gas.General.frcao(validData==1));
                handleDuration  = getLinFit(pointAboveCetCrit, validCEV_DS, cumsum(table.breathDuration(validData==1)));
            end
            %finding cumulative expired volume at CetCrit      
            CetCrit = table.criticalEndRatios(i);
            CEV = fzero(@(v) handleCet(v)-CetCrit,CEVmax/2);  % finding cumulative expired volume at CetCrit

            %Calculate the output variables
            FRC = handleFRC(CEV);
            LCI = CEV/FRC;
            duration=handleDuration(CEV);
            gas.ByFit.cet_norm(i)=CetCrit;
            gas.ByFit.lci(i)=LCI;
            gas.ByFit.frc(i)=FRC;
            gas.ByFit.nBreaths(i)=pointBelowCetCrit+length(find(validData(1:pointBelowCetCrit)==0));
            gas.ByFit.partialBreath(i)=(CEV-validCEV_DS(pointAboveCetCrit))/(validCEV_DS(pointBelowCetCrit)-validCEV_DS(pointAboveCetCrit));
            gas.ByFit.duration(i)=duration;
            gas.ByFit.handleCet{i}=handleCet;
            gas.ByFit.handleFRC{i}=handleFRC;
            gas.ByFit.minuteVentilation(i)=CEV/duration;
            gas.ByFit.breathRate(i)=gas.ByFit.nBreaths(i)/duration;
            gas.ByFit.deltaVolume(i)=validCEV_DS(pointBelowCetCrit)-CEV;
            
            if momentRatio
                 [m0, m1, m2] = getMomentRatios( ...
                                                validCet_norm(1:pointBelowCetCrit), ...
                                               [validCEV_DS(1:pointAboveCetCrit)./gas.General.frcao(1:pointAboveCetCrit); CEV/FRC]);
                gas.ByFit.momentRatio(i,1:3) = [m0, m1, m2];
            end
            
            data=gas.ByFit;
            fprintf('#breath=%3.1f: duration=%5.1f s, FRC=%6.1f ml, LCI=%5.2f at Cet(norm)=%4.1f %%, BR=%4.1f 1/min, MV=%6.1f ml/min\n', data.nBreaths(i)-(1-data.partialBreath(i)),data.duration(i),data.frc(i)*1e6,data.lci(i),CetCrit*100,data.breathRate(i)*60,data.minuteVentilation(i)*1e6*60);
        else
            fprintf('Cet(norm)=%4.1f %%: No valid breath found\n', CetCrit*100);
        end
    end
end

function handle = getLinFit(indexPointAbove, xData, yData)
           
        dx=xData(indexPointAbove+1)-xData(indexPointAbove);
        dy=yData(indexPointAbove+1)-yData(indexPointAbove);
        y0=mean([yData(indexPointAbove), yData(indexPointAbove+1)])-mean([xData(indexPointAbove), xData(indexPointAbove+1)])*dy/dx;
        handle  = @(v) y0+dy/dx*v;  
end

function handle = getExpFit(x, y)

    exp2 = @(a,x)(a(1).*exp(a(2)*x)+a(3).*exp(a(4)*x));     
    opts = optimset('MaxFunEvals',50000, 'MaxIter',10000);
    OLS = @(b) sum((exp2(b,x)-y).^2);
    minParams = fminsearch(OLS, [y(1); -0.5; y(1)/2; -0.5], opts);
    
    r2     = 1-sum((exp2(minParams,x)-y).^2);
    handle = @(x)(minParams(1)*exp(x*minParams(2))+minParams(3)*exp(x*minParams(4)));
end