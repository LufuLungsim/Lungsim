%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MATLAB script to test the BTPS correction
%
% Synopsis: two corresponding data files, a raw one and a BTPS corrected
% one by wBreath are specified compared to the own BTPS correction done
% within MATLAB
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 1.5, 19. Mai 2011
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fileName,pathName]=guiFileDialog(startDir,multi)

    if strcmp(multi,'off')
        windowText='Select data file';
    else
        windowText='Select one or more data files';
    end
    
    startDir=strcat(startDir,'\');
    [fileName,pathName] = uigetfile({'*.txt', 'wBreath data files (*.txt)'},windowText,startDir,'MultiSelect',multi);
end