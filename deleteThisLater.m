
% FilesStruct = {'Settings_Global'; 'Settings_Group'; 'Automatic'; 'Manual'};
% TabsStruct = {'Healty2'; 'Healthy3'; 'CF3'; 'CF3_2'; 'PCD'};
% 
% delayIVCO2   = cell(200,3);  
% delayCO2O2   = cell(200,3);
% delayCO2MMss = cell(200,3);
% FRC          = cell(200,3);
% LCI          = cell(200,3);
% 
% for j = 1:length(FilesStruct)
%     for i = 1:length(TabsStruct)        
%         rangeStart=1+(j-1)*50+(i-1)*10;
%         range=rangeStart:(rangeStart+9);
%         
%         delayIVCO2(range, 1)   = FilesStruct(j);
%         delayIVCO2(range, 2)   = TabsStruct(i);
%         delayIVCO2(range, 3)   = eval([FilesStruct{j} '.' TabsStruct{i} '(2:11,7)']);
%         
%         delayCO2O2(range, 1)   = FilesStruct(j); 
%         delayCO2O2(range, 2)   = TabsStruct(i);
%         delayCO2O2(range, 3)   = eval([FilesStruct{j} '.' TabsStruct{i} '(2:11,8)']);
%         
%         delayCO2MMss(range, 1) = FilesStruct(j); 
%         delayCO2MMss(range, 2) = TabsStruct(i);
%         delayCO2MMss(range, 3) = eval([FilesStruct{j} '.' TabsStruct{i} '(2:11,9)']);   
%         
%         FRC(range, 1)          = FilesStruct(j); 
%         FRC(range, 2)          = TabsStruct(i);
%         FRC(range, 3)          = eval([FilesStruct{j} '.' TabsStruct{i} '(2:11,30)']);
%         
%         LCI(range, 1)          = FilesStruct(j); 
%         LCI(range, 2)          = TabsStruct(i);
%         LCI(range, 3)          = eval([FilesStruct{j} '.' TabsStruct{i} '(2:11,31)']);
%     end
% end
% 
% % Files = {'Calibration_Settings_Global.xlsx';'Calibration_Settings_Per_Study_Group.xlsx';'Calibration_Automatic.xlsx';'Calibration_Manual.xlsx'};
% % FilesStruct = {'Settings_Global'; 'Settings_Group'; 'Automatic'; 'Manual'};
% % Tabs  = {'Healthy Set 2';'Healthy Set 3';'CF Set 3';'CF Set 3_2';'PCD'};
% % TabsStruct = {'Healty2'; 'Healthy3'; 'CF3'; 'CF3_2'; 'PCD'};
% % Path  = 'C:\Users\Jerry\Documents\MATLAB\Insel\Doku\Calibration Verification\';
% % 
% % for i = 1:length(Tabs)
% %     [~,~,A1] = xlsread([Path, Files{1}],Tabs{i}, 'A1:FG11');
% %     [~,~,A2] = xlsread([Path, Files{2}],Tabs{i}, 'A1:FG11');
% %     [~,~,A3] = xlsread([Path, Files{3}],Tabs{i}, 'A1:FG11');
% %     [~,~,A4] = xlsread([Path, Files{4}],Tabs{i}, 'A1:FG11');
% %      
% %     assignin('base', TabsStruct{i},     ...
% %               struct( FilesStruct{1}, {A1}, ...
% %                       FilesStruct{2}, {A2}, ...
% %                       FilesStruct{3}, {A3}, ...
% %                       FilesStruct{4}, {A4}));
% % end