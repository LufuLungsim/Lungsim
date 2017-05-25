%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MATLAB function to model temperature distribution within spirette
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 2.0, 30. August 2011
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [signal]=filterDelayTemp(signal,parameters)

    volumeBody2         =   parameters.Device.volumeBody2;
    volumeBody1         =   parameters.Device.volumeBody1;
    volumePrecap        =   parameters.Device.volumePrecap;
    volumePostcap       =   parameters.Device.volumePostcap;
    
    cumVolume           =   cumsum([0,volumeBody2,volumeBody1,volumePrecap,volumePostcap]);
    
    tauBody2            =   parameters.Device.tauBody2;
    tauBody1            =   parameters.Device.tauBody1;
    tauTemp             =   parameters.Device.tauTemp;
    
    tauArray            =   [tauBody1,tauBody2,tauTemp,tauTemp];
    
    tempBody2           =   parameters.Operation.Tbody2;
    tempBody1           =   parameters.Operation.Tbody1;
    tempSensor          =   parameters.Operation.Tsensor;
    
    tempBody            =   parameters.Operation.Tbody;
    tempBypass          =   parameters.Operation.TinsGas;
    
    tempArray           =   [tempBody1,tempBody2,tempSensor,tempSensor];
    
    dt                  =   parameters.Simulation.dt;
    
    t                   =   signal.ts;
    Iv                  =   signal.Iv;                      % flow
    
    nVtot               =	512;                            % # volume elements
    volumeTotal         =	cumVolume(end);                 % total modelled pre+post cap volume
    dVolume             =   volumeTotal/nVtot;              % volume increment
    volumeArray         =   [1:nVtot]*dVolume;              % array of cumulative volume elements
    
    partialVolume       =   cumVolume(end-1);
    iL                  =   floor(partialVolume/dVolume);
    weightL             =   partialVolume/dVolume-iL;
    weight              =   [weightL,1-weightL]';
    
    tauTemp = zeros(size(volumeArray));
    refTemp = zeros(size(volumeArray));
    for i = 2:length(cumVolume)
        tauTemp	= tauTemp+(volumeArray<=cumVolume(i)).*(volumeArray>cumVolume(i-1))*tauArray(i-1);
        refTemp	= refTemp+(volumeArray<=cumVolume(i)).*(volumeArray>cumVolume(i-1))*tempArray(i-1);
    end
    
    xi=1.0;
    Iv0=1e-4;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % shift type modelling including heat diffusion
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    inner=ones(length(t)+1,nVtot)*tempSensor;               % inner DOF of the delay line
    allNodes=[ones(length(t)+1,1)*tempBody,inner,ones(length(t)+1,1)*tempBypass];	% bc according 0 at left and 1 at right
    
    for n=2:length(Iv)+1                                    % temporal sweep
        delta=-(Iv(n-1)*dt)/dVolume;                        % scaled volume increment for n-th time step (convention: Iv<0==expiration)
        
        dj=floor(delta);                                    % full lattice spacing transport
        pj=delta-dj;                                        % partial lattice spacing transport
        j=2:1:nVtot+1;                                      % spatial sweep
        leftIndex =min(max(j-dj-1,1),nVtot+2);              % left spatial index including bc treatment
        rightIndex=min(max(j-dj,  1),nVtot+2);              % right spatial index including bc treatment
        
        allNodes(n,j)=(1-pj)*allNodes(n-1,rightIndex)+pj*allNodes(n-1,leftIndex);	% convection transport
        
        scalingIv=xi+(1+tanh(Iv(n-1)/Iv0))/2*(1-xi);
        diffFactorTemp=exp(-(dt)./(tauTemp*scalingIv));
        
        allNodes(n,2:end-1)=allNodes(n,2:end-1).*diffFactorTemp+refTemp.*(1-diffFactorTemp);    % heat diffusion to the walls
    end
    
    signal.delayTemp=allNodes(2:end,[1+iL,1+iL+1])*weight;	% interpolation of actual state
end
