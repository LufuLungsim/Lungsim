%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MATLAB function to model non-constant delays
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 1.0, 14. November 2014
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [delayLine]=singleStepDelay(delayLine,delta,nCompartments)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% shift type modelling
% with a routine that directly uses for loops (no MATLAB specific code)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deltaIndex=floor(delta);                                        % full lattice spacing transport
partialIndex=delta-deltaIndex;                              	% partial lattice spacing transport
newDof=delayLine;                                               % initialization of dof array
for j=2:nCompartments;                                          % delay line sweep
    leftIndex =min(max(j-deltaIndex-1,1),nCompartments);     	% left index including bc treatment
    rightIndex=min(max(j-deltaIndex,  1),nCompartments);      	% right index including bc treatment
    newDof(j)=(1-partialIndex)*delayLine(rightIndex)+partialIndex*delayLine(leftIndex);	% convective transport
end
delayLine=newDof;                                               % writing back to delayLine (might not be change in place!)

end

