%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculation Lagrange interpolation for signal fractional delays
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 1.0, 14. November 2014
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [b]=singleStepLagrange(delay,N)
%
% implements a lagrange interpolation filter of order N for fractional delay
%
% This implementation evaluates the filters for all time steps (delay can
% be a vector) in advance.
% For online-purposes, this filter coefficient routine should be called
% sequentially for each time step
%
Nfull=N*floor((delay)/N);       % separation of full orders of digital delays
Neff=Nfull+N;                   % effective filter order
b=zeros(Neff+1,1);            	% nominator coeffients
Deff=delay-Nfull;               % fractional delay
for k=0:N                       % calculation of lagrange interpolation coefficients
    nTerms=Deff-[0:N];
    nTerms(k+1)=1;
    dTerms=k-[0:N];
    dTerms(k+1)=1;
    b(Nfull+k+1)=prod(nTerms./dTerms);
end


