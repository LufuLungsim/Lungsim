%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optional smoothing of raw data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function signal = filterRawSignal(signal, nOrder, filterCutoff, dt) 
    if nOrder>0
        [a,b]           =   filterButterworth(filterCutoff,nOrder,dt); 	% discret butterworth filter by bilinear transformation
        signal.Iv       =   filterIIR(b,a,signal.Iv,[]);                % response to Iv signal
        signal.O2       =   filterIIR(b,a,signal.O2,[]);                % response to O2 signal
        signal.CO2      =   filterIIR(b,a,signal.CO2,[]);               % response to CO2 signal
        signal.MMss 	=   filterIIR(b,a,signal.MMss,[]);              % response to MMss signal
        signal.Ivss     =   filterIIR(b,a,signal.Ivss,[]);              % response to Iv signal
        signal.MMms     =   filterIIR(b,a,signal.MM,[]);                % response to MM signal
    end

    [a,b]               =   filterButterworth(1.2*2,2,dt);              % discret butterworth filter by bilinear transformation
    signal.IvssSmooth   =   filterIIR(b,a,signal.Ivss,[]);              % response to Iv signal
    firVector           =   ones(floor(0.6/dt),1);                     	% moving average mean
    firVector           =   firVector/length(firVector);
    signal.IvssSmooth1  =   filterFIR(firVector,signal.Ivss);
end