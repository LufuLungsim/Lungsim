%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Optimization of filter coefficients with non-linear Newton-Raphson
% algorithm
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 2.4, 20. Oktober 2011
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [coefficientsOptimal] = optimizeFilters(optimizingSet,coefficientScales,initialCoefficients,signal,parameters)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % setting simulation parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    relativeError	=   parameters.Simulation.relativeError*1.0;  	% relativ error iteration
    maximalIteration=   parameters.Simulation.maximalIteration;     % maximal iteration count
    verb            =   parameters.Simulation.verb;                 % verbosity level
    
    %method          =   'MMss+N2';                                	% optimization methods: MMss, N2, MMss+N2
    %method          =   'N2';                                      % optimization methods: MMss, N2, MMss+N2
    method          =   'MMss';                                     % optimization methods: MMss, N2, MMss+N2
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % optimization interation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nDOF            =   length(optimizingSet);
    delta           =   0.001;
    
    coefficients0	=   initialCoefficients;
    errorResiduals  =   1.0e99;
    errorIncrements =   1.0e99;
    nIteration      =   0;
    
    scaleN2         =   0.78;
    scaleMMss       =   0.03;
    
    errorResiduals0=1e99;
    errorIncrements0=1e99;
    while (nIteration<maximalIteration && errorResiduals>relativeError && errorIncrements>relativeError)
        dummy           =   signal;
        dummy           =   signalFilter(dummy,parameters,coefficients0);
        dummy          	=   signalCalc(dummy,parameters);
        switch(method)
            case 'MMss'
                f0              =  (dummy.MMssPoi3-dummy.MMssCalc);
            case 'N2'
                f0              =  (dummy.N2Poi3-mean(dummy.N2Poi3));
            case 'MMss+N2'
                f0              =   [(dummy.MMssPoi3-dummy.MMssCalc)/scaleMMss;(dummy.N2Poi3-mean(dummy.N2Poi3))/scaleN2];
            otherwise
                error(sprintf('ERROR: method %s unknown',method));
        end
        b               =  -f0;
        
        for i=1:nDOF
            coefficients1   =   coefficients0;
            coefficients1(optimizingSet(i))=coefficients1(optimizingSet(i))+coefficientScales(optimizingSet(i))*delta;
            dummy           =   signal;
            dummy           =   signalFilter(dummy,parameters,coefficients1);
            dummy          	=   signalCalc(dummy,parameters);
            switch(method)
                case 'MMss'
                    A(:,i)  		=   ((dummy.MMssPoi3-dummy.MMssCalc)-f0)/delta;
                case 'N2'
                    A(:,i)  		=   ((dummy.N2Poi3-mean(dummy.N2Poi3))-f0)/delta;
                case 'MMss+N2'
                    A(:,i)  		=   ([(dummy.MMssPoi3-dummy.MMssCalc)/scaleMMss;(dummy.N2Poi3-mean(dummy.N2Poi3))/scaleN2]-f0)/delta;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Linear solver based on the RQ decompostion (implicit in A\b)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        x   =   A\b;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Alternative linear solver that aviods the RQ decompostion implicit in A\b
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     Ab  =   A'*b;
        %     AtA =   A'*A;
        %     x   =   inv(AtA)*Ab;
        
        deltaCoefficients   =   zeros(size(coefficients0));
        for i=1:nDOF
            deltaCoefficients(optimizingSet(i))=coefficientScales(optimizingSet(i))*x(i);
        end
        
        errorIncrements =   norm(deltaCoefficients)/length(coefficients0);
        switch(method)
            case 'MMss'
                errorResiduals  =   norm(b)/norm(dummy.MMssPoi3);
            case 'N2'
                errorResiduals  =   norm(b)/norm(dummy.N2Poi3);
            case 'MMss+N2'
                errorResiduals  =   norm(b)/sqrt(norm(dummy.MMssPoi3)^2+norm(dummy.N2Poi3)^2);
        end
        
        if errorResiduals>errorResiduals0
            fprintf('stopping iteration due to error increase at: ErrorResiduals=%e, ErrorIncrements=%e, #iteration=%d\n',errorResiduals0,errorIncrements0,nIteration);
            break;
        end
        errorResiduals0=errorResiduals;
        errorIncrements0=errorIncrements;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Updating coefficients
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        coefficients0	=   coefficients0+deltaCoefficients*1.00;	% possiblility for overrelaxation to improve stablity, but proves not better
        coefficients0(3)=   max(coefficients0(3),0.0000001);         % remove negative filter coefficient
        
        if verb
            fprintf('%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t\n',coefficients0);
        end
        
        fprintf('ErrorResiduals=%e, ErrorIncrements=%e, #iteration=%d\n',errorResiduals,errorIncrements,nIteration);
        
        nIteration      =   nIteration+1;
    end
    
    if nIteration==maximalIteration
        warning('maximal iteration count reached: possible bad convergence');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % returning optimal coefficients
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    coefficientsOptimal =   coefficients0;
end