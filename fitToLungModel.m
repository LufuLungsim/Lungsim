function gas = fitToLungModel(table, gas, parameters)
    addpath('.\FMINSEARCHBND');
    graphState = parameters.Simulation.graphState;
    
    %Prepare the used variables
    vExp    = table.volExp;
    vInsp   = table.volInsp;
    vReInsp = gas.reInsp;
    vInsp   = vInsp + vReInsp;
    cExp    = gas.General.cet_norm;
    cInsp   = vReInsp;
    c0      = gas.General.cet_start(1);
    tInsp   = table.inspDuration;
    tExp    = table.expDuration;
    tBreath = tInsp+tExp;
    FRC     = mean(gas.General.frcao(end-3:end));
    meanTidal = mean(table.volTidal);
    vLung   = FRC + mean(table.volTidal);
    opts    = optimset('MaxFunEvals',5000000, 'MaxIter',1000000); %Maybee add TolFun TolX
 
    %Set starting and boundary values to fit the two compartment model
    %The values are the lung volume and the proportion of the volume of the
    %turbulent compartment as well as the Diffusion barrier factor. The
    %volume of the turbulent and diffusion compartment are calculated based
    %on the lung volume and proportion factor
    x0 = [vLung; 0.5; 0.001];
    LB = [vLung-meanTidal; 0.1; 1e-9];
    UB = [vLung+0.5*meanTidal; 0.9; inf];
    
    %Fit the model to the measured Data
    OLS = @(b) testTwoCompartmentModel(b(1), b(2), b(3), vExp, vInsp, c0, cExp, cInsp, tBreath);
    twoCompartments   = fminsearchbnd(OLS,x0,LB,UB,opts);
        
    %Use the found parameters to evaluate the function after some math
    vTurbulent = twoCompartments(1)*twoCompartments(2);
    vDiffusive = twoCompartments(1)-vTurbulent;
    cetFit     = fitTwoCompartmentModel(vTurbulent, vDiffusive, twoCompartments(3), vExp, vInsp, c0, cExp, cInsp, tBreath);
    r2         = getR2y(cExp, cetFit);
    gas.twoCompartments = struct('vTurbulent', vTurbulent, 'vDiffusion', vDiffusive, 'D', twoCompartments(3), 'r2', r2);
    
    if graphState == 1
        getOrMakeFigure('Conduction and diffusion compartment Approximation 2 Chamber Model');
        hold off
        plot(cExp, '-r')
        hold on
        plot(cetFit, '-b')
        legend('Measured Values', 'Approximation by Model');
        xlabel('Breath Nr.'); ylabel('Gas concentration');
    end
    
    
    %Set starting and boundary values to fit the two compartment model
    %The values are the lung volume and the proportion of the volume of the
    %turbulent compartment as well as the Diffusion barrier factor. The
    %volume of the turbulent and diffusion compartment are calculated based
    %on the lung volume and proportion factor
    x0 = [vLung; 1/3; 0.1; 0.001];
    LB = [vLung-0.5*meanTidal; min(1.5*meanTidal/vLung, 0.5); 0.09; 1e-9];
    UB = [vLung+0.5*meanTidal; 0.9                          ; 0.11; inf];
    
    %Fit the model to the measured Data
    OLS = @(b) testThreeCompartmentModel(b(1), b(2), b(3), b(4),  vExp, vInsp, c0, cExp, cInsp, tBreath);
    threeCompartments   = fminsearchbnd(OLS,x0,LB,UB,opts);
    
    %Use the found parameters to evaluate the function after some math
    vTurbulent = threeCompartments(1) * threeCompartments(2);
    vLaminar   = threeCompartments(1) * threeCompartments(3);
    vDiffusive = threeCompartments(1) - vTurbulent - vLaminar;
    cetFit     = fitThreeCompartmentModel(vTurbulent, vLaminar, vDiffusive, threeCompartments(4), vExp, vInsp, c0, cExp, cInsp, tBreath);
    r2         = getR2y(cExp, cetFit);
    gas.threeCompartments = struct('vTurbulent', vTurbulent, 'vLaminar', vLaminar, 'vDiffusion', vDiffusive, 'D', threeCompartments(4), 'r2', r2);
    
    if graphState == 1
        getOrMakeFigure('Conduction and diffusion compartment Approximation 3 Chamber Model');
        hold off
        plot(cExp, 'r-')
        hold on
        plot(cetFit, 'b-')
        legend('Measured Values', 'Approximation by Model');
        xlabel('Breath Nr.'); ylabel('Gas concentration');
    end
    
end

%Test the two compartment model
%
%Could also be done directly by the OSL function but its a little easier to
%read this way
function fitMeasure = testTwoCompartmentModel( vLung, pTurbulent, D, vExp, vInsp, c0, cExp, cInsp, tBreath)
    vTurbulent=vLung*pTurbulent;
    vDiffusive=vLung-vTurbulent;
    y = fitTwoCompartmentModel(vTurbulent, vDiffusive, D, vExp, vInsp, c0, cExp, cInsp, tBreath);
    
    if isreal(y)
        curveFit   = sum(((y-cExp)./cExp).^2);
    else
        curveFit   = 1e9;
    end
    fitMeasure = curveFit;
end

%Test the three compartment model
%
%
%Additionally to the fit of the data a check is made to ensure that the sum
%of the compartment volume is approximately the same as the total lung
%volume
function fitMeasure = testThreeCompartmentModel( vLung, pTurbulent, pLaminar, D, vExp, vInsp, c0, cExp, cInsp, tBreath)
    vTurbulent = vLung * pTurbulent;
    vLaminar   = vLung * pLaminar;
    vDiffusive = vLung - vTurbulent - vLaminar;
    
    if vDiffusive/vLung > 0.1
        y = fitThreeCompartmentModel(vTurbulent, vLaminar, vDiffusive, D, vExp, vInsp, c0, cExp, cInsp, tBreath);
    
        curveFit   = sum(((y-cExp)./cExp).^2);
        volumeFit  = 0;%(pTurbulent+pLaminar+vDiffusive/vLung-1)^2;
    
        fitMeasure = curveFit+volumeFit;
    else
        fitMeasure = 1e9;
    end
end

%
% Calculate the diffusion of the gas concentration based on a simple model
% of the lung consisting of two compartments. The first compartment which
% is connected to the outside world is ventilated by breathing, assuming a 
% perfect mixture of the incoming gas with the gas mixture in the compartment.  
% A second compartment, connected to the first but seperated by a membrane
% which is permable by gas and shows diffusive properties holds a reservoir
% of gas. 
%
%                       Lungmodel
%                 _________________________
%             |__|            ¦            |
%  outside    ___   Turbulent ¦  Diffusive |
%             |  |____________¦____________| 
%
function cExpCalc = fitTwoCompartmentModel(vTurbulent, vDiffusive, D, vExp, vInsp, c0, cExp, cInsp, tBreath)
   
    b=length(cExp);
                                                                                                             
    cDiffusive    = zeros(b+1,1);                                                 
    cTurbulent    = zeros(b+1,1);                                                
    cDiffusive(1) = c0;
    cTurbulent(1) = c0;
    
    for i=2:b+1    
        cDiffusive(i) = cDiffusive(i-1)*(1+tBreath(i-1)*D*log(cTurbulent(i-1)/cDiffusive(i-1)));                
        cTurbulent(i)   = ((cTurbulent(i-1)*(vTurbulent-vExp(i-1))+cInsp(i-1)*vInsp(i-1))...
                    +(cDiffusive(i-1)-cDiffusive(i))*vDiffusive)/vTurbulent;
    end
    
    text(1, 0.05, ['vTurbulent: ' num2str(vTurbulent) ' vDiffusive: ' num2str(vDiffusive) ' ']);
    cExpCalc = cTurbulent(2:end);
end


%
% Calculate the diffusion of the gas concentration based on a simple model
% of the lung consisting of three compartments. The first compartment which
% is connected to the outside world is ventilated by breathing, assuming a 
% perfect mixture of the incoming gas with the gas mixture in the compartment.
% A second compartment between the Turbulent and Diffusive compartment as
% seen in the two compartment method is added. This compartment is assumed
% to have laminar flow and always has the average tracer gas concentration
% of the first and third compartment. 
% A third compartment, connected to the second but seperated by a membrane
% which is permable by gas and shows diffusive properties, holds a reservoir
% of gas. 
%
%                              Lungmodel
%                 _________________________________
%                |           _________¦            |
%                | Turbulent _Laminar_¦  Diffusive |
%             |__|           _________¦            |
%  outside     __            _________¦            |
%             |  |____________________¦____________| 
%
function cExpCalc = fitThreeCompartmentModel(vTurbulent, vLaminar, vDiffusion, D, vExp, vInsp, c0, cExp, cInsp, tBreath)
   
    b=length(cExp);
                                                                                                             
    cDiffusion = zeros(b+1,1);
    cLaminar   = zeros(b+1,1);
    cTurbulent = zeros(b+1,1); 
    
    mDiffusion = zeros(b+1,1);
    mLaminar   = zeros(b+1,1);
    mTurbulent = zeros(b+1,1);
    
    cDiffusion(1) = c0;
    cLaminar(1)   = c0;
    cTurbulent(1) = c0;
    
    mDiffusion(1) = cDiffusion(1)*vDiffusion;
    mLaminar(1)   = cLaminar(1)*vLaminar;
    mTurbulent(1) = cTurbulent(1)*vTurbulent;

    for i=2:b+1   
        cDiffusion(i) =  cDiffusion(i-1)*(1+tBreath(i-1)*D*log(cLaminar(i-1)/cDiffusion(i-1))); 
        %mDiffusion(i) =  cDiffusion(i-1)*(1+tBreath(i-1)*D*log(cLaminar(i-1)/cDiffusion(i-1)))*vDiffusion;
        mDiffusion(i) =  cDiffusion(i)*vDiffusion;
        
        cLaminar(i)   =  (cDiffusion(i)+cTurbulent(i-1))/2;
        mLaminar(i)   =  cLaminar(i)*vLaminar;
        
        
        
        cTurbulent(i)   = ((cTurbulent(i-1)*(vTurbulent-vExp(i-1))+cInsp(i-1)*vInsp(i-1))...
                           +(cLaminar(i-1)-cLaminar(i))*vLaminar)/vTurbulent;
%         mTurbulent(i) =  cTurbulent(i-1)*(vTurbulent-vExp(i-1))  ...%Amount of gas in the compartment in the last breath minus the expired amount
%                          +cInsp(i-1)*vInsp(i-1)                  ...%Plus the reinspired gas
%                          -(mLaminar(i)-mLaminar(i-1));           ...%Plus the amount going to/coming from the laminar compartment   
        mTurbulent(i) = cTurbulent(i)*vTurbulent;
        
    end
    
    text(1, 0.05, ['vTurbulent: ' num2str(vTurbulent) ' vLaminar: ' num2str(vLaminar) ' vDiffusive: ' num2str(vDiffusion) ' ']);
    cExpCalc = cTurbulent(2:end);
end

function r2 = getR2y(y, yf)
    r2     = 1-sum((y-yf).^2)/sum((y-mean(y)).^2);
end