%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MATLAB function to detect breaths
% Attention: This implementation relays on the delayLine signal. This
% signal was introduced to enable temperature- and humidity induced error
% correction for the MMms-signal. It proves, however, that the delayLine
% changed between 0->1 and 1->0 almost ideally mark the change from
% inspiration to expiration and vice-versa.
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 4.0, 30. August 2011
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%returns index (index of the start of breaths in flow) and breaths (timestamp of the breaths as in time) 
function [index,breaths]=breathDetection(time,flow,parameters)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % setup of delay line with parameters different from spirette settings
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dt                  =   parameters.Simulation.dt;
    nVtot               =	256;                            % # volume elements
    volumeTotal         =	10e-6;                          % pre+post cap volume
    if parameters.Simulation.MissPlexy == 1 || parameters.Simulation.Set1 == 1
       minBreathVolume  =   50e-7; 
    else 
        minBreathVolume =   50e-6;                          % minimal breath volume (+ or -)
    end
    dVolume             =   volumeTotal/nVtot;              % volume increment
    iL                  =   floor(volumeTotal/2/dVolume);
    weightL             =   volumeTotal/2/dVolume-iL;
    weight              =   [weightL,1-weightL]';
    
    inner       =   ones(length(time)+1,nVtot);             % inner DOF of the delay line
    allNodes=[zeros(length(time)+1,1),inner,ones(length(time)+1,1)];	% bc according 0 at left and 1 at right
    for n=2:length(flow)+1                                  % temporal sweep
        delta=-(flow(n-1)*dt)/dVolume;                      % scaled volume increment for n-th time step (convention: Iv<0==expiration)
        
        dj=floor(delta);                                    % full lattice spacing transport
        pj=delta-dj;                                        % partial lattice spacing transport
        j=2:1:nVtot+1;                                      % spatial sweep
        leftIndex =min(max(j-dj-1,1),nVtot+2);              % left spatial index including bc treatment
        rightIndex=min(max(j-dj,  1),nVtot+2);              % right spatial index including bc treatment
        allNodes(n,j)=(1-pj)*allNodes(n-1,rightIndex)+pj*allNodes(n-1,leftIndex);	% convection transport
    end
    delay	=   allNodes(2:end,[1+iL,1+iL+1])*weight;       % interpolation of actual state
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % finding all indices in delay with a sign change
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    i=0;
    j=1;
    indexRaw=[];
    actualSign=sign(delay(1)-0.5);
    while j<length(delay)
        newSign=sign(delay(j)-0.5);
        if ~(newSign==actualSign || newSign==0)
            i=i+1;
            indexRaw(i)=j;
            actualSign=newSign;
        end
        j=j+1;
    end
    indexRaw=[indexRaw,length(delay)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % adjusting indices for actual flow sign change (to remove effects of delay
    % line filtering)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    index=indexRaw;
    for i=1:length(indexRaw)
        j=indexRaw(i)-1;
        while(j>0)
            if(sign(flow(j+1))~=sign(flow(j)))
                break;
            end
            j=j-1;
        end
        index(i)=j;
    end
    if index(1)==0                              % remove double indices at start and end
        index=index(2:end);
    end
    if index(end)==index(end-1)
        index=index(1:end-1);  
    end                        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % remove breaths with too small volumes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j=1;
    newIndex(1)=index(1);
    i=2;
    while i<length(index)       % find initial breath
        subIndex=(index(i-1):1:index(i));
        vol =trapz(time([subIndex]),flow([subIndex]));
        if abs(vol)>minBreathVolume
            break;
        end
        i=i+1;
    end
    while i<length(index)       % iteration
        subIndex=(index(i-1):1:index(i));
        vol = trapz(time([subIndex]),flow([subIndex]));
        if abs(vol)>minBreathVolume
            j=j+1;
            newIndex(j)=index(i);
            i=i+1;
        else
            newIndex(j)=index(i+1);
            i=i+2;
        end
    end
    
    index=newIndex;
    
    breaths = [];                               % determine breath times
    for i=1:length(index)-1
        breaths(i)=time(index(i+1));
    end
    
end


