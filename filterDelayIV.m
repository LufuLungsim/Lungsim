%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MATLAB function to model non-constant delays
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 2.0, 30. August 2011
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [signal]=filterDelayIV(signal,parameters)

    volumePrecap        =   parameters.Device.volumePrecap;
    volumePostcap       =   parameters.Device.volumePostcap;
    dt                  =   parameters.Simulation.dt;
    
    t                   =   signal.ts;
    Iv                  =   signal.Iv;                      % flow
    
    nVtot               =	256;                             % # volume elements
    volumeTotal         =	volumePrecap + volumePostcap;	% pre+post cap volume
    dVolume             =   volumeTotal/nVtot;              % volume increment
    
    iL                  =   floor(volumePrecap/dVolume);
    weightL             =   volumePrecap/dVolume-iL;
    weight              =   [weightL,1-weightL]';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % shift type modelling
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    inner       =   ones(length(t)+1,nVtot);                % inner DOF of the delay line
    allNodes=[zeros(length(t)+1,1),inner,ones(length(t)+1,1)];	% bc according 0 at left and 1 at right
    
    for n=2:length(Iv)+1                                    % temporal sweep
        delta=-(Iv(n-1)*dt)/dVolume;                        % scaled volume increment for n-th time step (convention: Iv<0==expiration)
        
        dj=floor(delta);                                    % full lattice spacing transport
        pj=delta-dj;                                        % partial lattice spacing transport
        j=2:1:nVtot+1;                                      % spatial sweep
        leftIndex =min(max(j-dj-1,1),nVtot+2);              % left spatial index including bc treatment
        rightIndex=min(max(j-dj,  1),nVtot+2);              % right spatial index including bc treatment
        allNodes(n,j)=(1-pj)*allNodes(n-1,rightIndex)+pj*allNodes(n-1,leftIndex);	% convection transport
    end
    
    signal.delayLine    =   allNodes(2:end,[1+iL,1+iL+1])*weight;	% interpolation of actual state
end    
