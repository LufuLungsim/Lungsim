function varargout = guiWbreath(varargin)
% GUIWBREATH M-file for guiWbreath.fig
%      GUIWBREATH, by itself, creates a new GUIWBREATH or raises the existing
%      singleton*.
%
%      H = GUIWBREATH returns the handle to a new GUIWBREATH or the handle to
%      the existing singleton*.
%
%      GUIWBREATH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIWBREATH.M with the given input arguments.
%
%      GUIWBREATH('Property','Value',...) creates a new GUIWBREATH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUIWBREATH before guiMassSpec_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to guiWbreath_OpeningFcn via varargin.
%
%      *See GUIWBREATH Options on GUIDE's Tools menu.  Choose "GUIWBREATH allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help guiWbreath

% Last Modified by GUIDE v2.5 03-Jul-2015 09:27:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @guiWbreath_OpeningFcn, ...
                   'gui_OutputFcn',  @guiWbreath_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% Using non-native system dialogs as workaround for an error occuring
% for large directory listings in uigetfile....
setappdata(0,'UseNativeSystemDialogs',false);

% --- Executes just before guiWbreath is made visible.
function guiWbreath_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to guiWbreath (see VARARGIN)

% Choose default command line output for guiWbreath
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes guiWbreath wait for user response (see UIRESUME)
% uiwait(handles.figure1);

if isempty(get(hObject,'UserData'))
    parameters=setParametersWbreath('parametersWbreath.mat',1);

    defaultParameters=setParametersWbreath('parametersWbreath.mat',0);
    oldVersionString=parameters.Simulation.versionString;
    newVersionString=defaultParameters.Simulation.versionString;
    
    foundPattern=findstr(oldVersionString,'.');
    if length(foundPattern)>1
        oldVersionString=oldVersionString(1:foundPattern(2)-1);
    end
    foundPattern=findstr(newVersionString,'.');
    if length(foundPattern)>1
        newVersionString=newVersionString(1:foundPattern(2)-1);
    end
    
    if ~strcmp(oldVersionString,newVersionString)
        warningString=sprintf('%s of parameters.mat is older than %s! Default settings used instead!',oldVersionString,newVersionString);
        handle=warndlg(warningString,'Warning (parameters.mat)','non-modal');
        uiwait(handle);
        parameters=defaultParameters;
    else
        parameters.Simulation.versionString=defaultParameters.Simulation.versionString;
    end
    
    set(handles.output,'UserData',parameters);
end
initializeParameters(hObject, handles)



% --- Outputs from this function are returned to the command line.
function varargout = guiWbreath_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on selection change in fileList.
function fileList_Callback(hObject, eventdata, handles)
% hObject    handle to fileList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns fileList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fileList
try
    interfaceObject=findobj(handles.output,'Enable','on'); % disable interface during evaluation
    set(interfaceObject,'Enable','off');
    drawnow;
    
    sourcePathName=get(handles.filePath,'String');
    choice=get(hObject,'Value');
    fileList=get(hObject,'String');
    if ~iscell(fileList)
        fileList={fileList};
    end
    fileName=fileList(choice);
    parameters=get(handles.output,'UserData');
    
    destinationPathName=get(handles.resultDirPath,'String');
    if ~exist(destinationPathName,'dir')
        buttonString = questdlg('Destination path does no longer exist...!','Warning','Choose new path','Use default','Use default');
        if strcmp(buttonString,'Use default')
            destinationPathName=getenv('TEMP');
        else
            destinationPathName=uigetdir(sourcePathName,'select result directory'); % starts at the present choice file directory!
            if isequal(destinationPathName,0)
                warning('(fileList_Callback): setting default destination path: no regular user choice!')
                destinationPathName=getenv('TEMP');
            end
        end
        set(handles.resultDirPath,'String',destinationPathName)
        parameters.Simulation.resultDirPath=destinationPathName;
    end
    
    parameters=applyCorrectionWbreath(sourcePathName,fileName,destinationPathName,parameters,1);
    
    set(handles.output,         'UserData',parameters);
    
catch err
    fclose('all');
    warning('(fileList_Callback): %s (in %s on line %g)\n', err.message,err.stack(1).name,err.stack(1).line);
end % try/catch

set(interfaceObject,'Enable','on');



% --- Executes on button press in choiceButton.
function choiceButton_Callback(hObject, eventdata, handles)
% hObject    handle to choiceButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%showSignals;
%cleanup(0); % must improve choice of figures here
try
    interfaceObject=findobj(handles.output,'Enable','on'); % disable interface during evaluation
    set(interfaceObject,'Enable','off');
    drawnow;
    
    parameters=get(handles.output,'UserData');
    
    pathName=get(handles.filePath,'String');
    if ~exist(pathName,'dir')
        if exist(parameters.Simulation.directoryDatafiles,'dir')
            pathName=parameters.Simulation.directoryDatafiles;
        else
            pathName=pwd;
        end
    end
    [fileName,pathName] = guiFileDialog(pathName,'on');
    if ~isequal(pathName,0)
        parameters.Simulation.directoryDatafiles=pathName;
    end
    if isequal(fileName,0)
        set(interfaceObject,'Enable','on');
        return;
    end
    set(handles.filePath,'String',pathName)
    set(handles.fileList,'String',fileName)
    set(handles.fileList,'Value',1)             % reset choice to 1st entry.
    
    set(handles.output,'UserData',parameters);

catch err
    fclose('all');
    warning('(choiceButton_Callback): %s (in %s on line %g)\n', err.message,err.stack(1).name,err.stack(1).line);
end % try/catch

set(interfaceObject,'Enable','on');



% --- Executes during object creation, after setting all properties.
function fileList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fileList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close all;
clear;



% --- Executes on button press in verbosityButton.
function verbosityButton_Callback(hObject, eventdata, handles)
% hObject    handle to verbosityButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of verbosityButton
verbosityState=get(hObject,'Value');
parameters=get(handles.output,'UserData');
parameters.Simulation.verb=verbosityState;
set(handles.output,'UserData',parameters);
fprintf('verbosity is equal %d\n',verbosityState)



% --- Executes on button press in quitButton.
function quitButton_Callback(hObject, eventdata, handles)
% hObject    handle to quitButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
parameters=get(handles.output,'UserData');
parameters.Calibration.settingsFilePath=strcat(pwd,'\','parametersWbreath.mat');
save('parametersWbreath.mat','parameters');
closereq



% --- Executes on button press in closeButton.
function closeButton_Callback(hObject, eventdata, handles)
% hObject    handle to closeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% closing all figures (including MATLAB type figures)
% should be improved to include only children
figureHandles = findall(0,'type','figure');
for h=figureHandles'
    if isempty(strfind(get(h,'name'),'LungSim for '))
        close(h)
    end
end



% --- Executes on button press in loadParameter.
function loadParameter_Callback(hObject, eventdata, handles)
% hObject    handle to loadParameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
parameters=get(handles.output,'UserData');

fullPathName=get(handles.settingsPath,'String');
if ~exist(fullPathName,'dir')
    pathName=strcat(pwd,'\*.mat');
end
windowText='Select calibration parameter file';
startDir=fullPathName;
[fileName,pathName] = uigetfile({'*.mat', 'calibration files (*.mat)'},windowText,startDir,'MultiSelect','off');
if ~isequal(fileName,0)
    fullName=strcat(pathName,fileName);
    if exist(fullName,'file')
        actualParameters=parameters;
        load(fullName);
        
        loadedVersionString=parameters.Simulation.versionString;
        actualVersionString=actualParameters.Simulation.versionString;
    
        if ~strcmp(loadedVersionString,actualVersionString)
            warningString=sprintf('%s of parameters.mat is older than %s!\nMissing defaults are updated!\n',loadedVersionString,actualVersionString);
            handle=warndlg(warningString,'Warning (parameters.mat)','non-modal');
            uiwait(handle);
            
            names1=fieldnames(actualParameters);
            for i=1:length(names1)
                names2=fieldnames(actualParameters.(names1{i}));
                for j=1:length(names2)
                    if ~isfield(parameters.(names1{i}),names2{j})
                        fprintf('(loadParameter) updating %s.%s=%s\n',names1{i},names2{j},mat2str(actualParameters.(names1{i}).(names2{j})));
                        parameters.(names1{i}).(names2{j})=actualParameters.(names1{i}).(names2{j});
                    end
                end
            end
            parameters.Simulation.versionString=actualParameters.Simulation.versionString;	% update of version string
        end
    end
    set(handles.settingsPath,'String',parameters.Calibration.settingsFilePath)
else
    fprintf('cannot load calibration parameter file\n')
end

set(handles.output,'UserData',parameters);

initializeParameters(hObject, handles)



% --- Executes on button press in saveButton.
function saveButton_Callback(hObject, eventdata, handles)
% hObject    handle to saveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
parameters=get(handles.output,'UserData');

fileName=strcat('calibration_',date);
fullName=strcat(pwd,'\',fileName);
[fileName,pathName] = uiputfile({'*.mat', 'calibration files (*.mat)'},'save settings file name',fullName);
if isequal(fileName,0)
    fprintf('cannot save calibration parameter file, no file chosen\n')    
else
    fullName=strcat(pathName,fileName);
    parameters.Calibration.settingsFilePath=fullName;
    save(fullName,'parameters');
    set(handles.settingsPath,'String',fullName);
    set(handles.output,'UserData',parameters);
end



function initializeParameters(hObject, handles)
% initializes parameters from file
parameters=get(handles.output,'UserData');

backgroundColor=get(0,'defaultUicontrolBackgroundColor');
foregroundColor=get(0,'defaultUicontrolForegroundColor');

set(handles.vPre,                       'String',num2str(parameters.Device.volumePrecap))
set(handles.vPost,                      'String',num2str(parameters.Device.volumePostcap))

set(handles.lowerBoundEE,               'String',num2str(parameters.Simulation.lowerBoundEE*100))
set(handles.upperBoundEE,               'String',num2str(parameters.Simulation.upperBoundEE*100))
set(handles.lowerBoundEI,               'String',num2str(parameters.Simulation.lowerBoundEI*100))
set(handles.upperBoundEI,               'String',num2str(parameters.Simulation.upperBoundEI*100))

set(handles.nBreathMax,                 'String',num2str(parameters.Simulation.nBreathMax))

set(handles.verbosityButton,            'Value', parameters.Simulation.verb)
set(handles.onlyPlotButton,             'Value', parameters.Simulation.onlyPlotState)
set(handles.graphButton,                'Value', parameters.Simulation.graphState)

set(handles.SF6Button,                  'Value', parameters.Simulation.SF6)
set(handles.HeButton,                   'Value', ~parameters.Simulation.SF6)
set(handles.sideStreamButton,           'Value', parameters.Simulation.sideStream)
if parameters.Simulation.SF6 == 1
    set(handles.sideStreamButton,'Enable','off');
else
    set(handles.sideStreamButton,'Enable','on');
end

set(handles.nBreathTidal,               'String',num2str(parameters.Simulation.nBreathTidal))
set(handles.fromBreathTidal,            'String',num2str(parameters.Simulation.fromBreathTidal))
set(handles.tidalMeanButton,            'Value', parameters.Simulation.tidalMeanState)

set(handles.tidalButton,                'Value', parameters.Simulation.tidalState)
set(handles.croppingButton,             'Value', parameters.Simulation.croppingState)
set(handles.washoutButton,              'Value', parameters.Simulation.washout)
set(handles.washinButton,               'Value', ~parameters.Simulation.washout)
if parameters.Simulation.croppingState || parameters.Simulation.tidalState
    set(handles.washoutButton,          'Enable', 'Off')
    set(handles.washinButton,           'Enable', 'Off')
else
    set(handles.washoutButton,          'Enable', 'On')
    set(handles.washinButton,           'Enable', 'On')    
end
if parameters.Simulation.tidalState
%     set(handles.minVentButton,          'Enable', 'Off')
    set(handles.lciStandardButton,      'Enable', 'Off')
    set(handles.lciByFitButton,         'Enable', 'Off')
    set(handles.TOevaluationButton,     'Enable', 'Off')
%    set(handles.momentRatiosButton,     'Enable', 'Off')
else
%     set(handles.minVentButton,          'Enable', 'On')
    set(handles.lciStandardButton,      'Enable', 'On')
    set(handles.lciByFitButton,         'Enable', 'On')
    set(handles.TOevaluationButton,     'Enable', 'On')
%    set(handles.momentRatiosButton,     'Enable', 'Off')
end

set(handles.meanButtonEE,               'Value', parameters.Simulation.meanEE==1)
set(handles.medianButtonEE,             'Value', parameters.Simulation.meanEE==2)
set(handles.byFitButtonEE,              'Value', parameters.Simulation.meanEE==3)
set(handles.meanButtonEI,               'Value', parameters.Simulation.meanEI==1)
set(handles.medianButtonEI,             'Value', parameters.Simulation.meanEI==2)
set(handles.byFitButtonEI,              'Value', parameters.Simulation.meanEI==3)
set(handles.rmsCrit,                    'String',num2str(parameters.Simulation.rmsCrit*100)) % ratio as % value

set(handles.fixedDeltaButton,           'Value', parameters.Simulation.deltaMMState)
set(handles.deltaMMFix,                 'String',num2str(parameters.Simulation.deltaMMFixed*1000)) % ratio as % value
if parameters.Simulation.deltaMMState
	set(handles.deltaMMFix,             'Enable','On')
else
	set(handles.deltaMMFix,             'Enable','Off')
end

set(handles.rmsCrit,                    'String',num2str(parameters.Simulation.rmsCrit*100)) % ratio as % value
set(handles.interactiveYesButton,       'Value', parameters.Simulation.interactiveYes)
set(handles.interactiveNoButton,        'Value', parameters.Simulation.interactiveNo)
set(handles.interactiveDeviationButton,	'Value', parameters.Simulation.interactiveDeviation)

if parameters.Simulation.interactiveNo 
    set(handles.rmsCrit,'Background',backgroundColor,'Foreground',[0.5 0.5 0.5]);
else
    set(handles.rmsCrit,'Background','white','Foreground',foregroundColor);
end

set(handles.minVentButton,              'Value', parameters.Simulation.minVentState)
set(handles.lciStandardButton,      	'Value', parameters.Simulation.lciStandard)
set(handles.lciByFitButton,          	'Value', parameters.Simulation.lciByFit)
set(handles.TOevaluationButton,      	'Value', parameters.Simulation.TOevaluation)
set(handles.momentRatiosButton,      	'Value', parameters.Simulation.momentRatios)

if ~strcmp('',parameters.Calibration.settingsFilePath)
    set(handles.settingsPath,	        'String',parameters.Calibration.settingsFilePath);
end

set(handles.versionString,              'String',parameters.Simulation.versionString)
set(handles.resultDirPath,              'String',parameters.Simulation.resultDirPath);

set(handles.output,'UserData',parameters);



% --- Executes on button press in resultDirButton.
function resultDirButton_Callback(hObject, eventdata, handles)
% hObject    handle to resultDirButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
parameters=get(handles.output,'UserData');

fullName=strcat(pwd);
pathName=uigetdir(fullName,'select result directory');
if ~isequal(pathName,0)
    set(handles.resultDirPath,'String',pathName)
    parameters.Simulation.resultDirPath=pathName;
end

set(handles.output,'UserData',parameters);



% --- Executes on button press in batchProcessingButton.
function batchProcessingButton_Callback(hObject, eventdata, handles)
% hObject    handle to batchProcessingButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    interfaceObject=findobj(handles.output,'Enable','on'); % disable interface during evaluation
    set(interfaceObject,'Enable','off');
    drawnow;
    
    parameters=get(handles.output,'UserData');
    
    sourcePathName=get(handles.filePath,'String');
    fileList=get(handles.fileList,'String');
    if ~iscell(fileList)
        fileList={fileList};
    end
    
    if strcmp('Signal Files to Process',cell2mat(fileList(1)))
        fprintf('INFO: no files chosen / present; use <Data File Choice> first\n')
        return;
    end
    
    destinationPathName=get(handles.resultDirPath,'String');
    if ~exist(destinationPathName,'dir')
        buttonString = questdlg('Destination path does no longer exist...!','Warning','Choose new path','Use default','Use default')
        if strcmp(buttonString,'Use default')
            destinationPathName=getenv('TEMP');
        else
            destinationPathName=uigetdir(sourcePathName,'select result directory'); % starts at the present choice file directory!
            if isequal(destinationPathName,0)
                warning('(fileList_Callback): setting default destination path: no regular user choice!')
                destinationPathName=getenv('TEMP');
            end
        end
        set(handles.resultDirPath,'String',destinationPathName)
        parameters.Simulation.resultDirPath=destinationPathName;
        set(handles.output,'UserData',parameters);
    end
    
    graphState=get(handles.graphButton,'Value');                                        % read user settings graph state
    parameters.Simulation.verb=parameters.Simulation.verb*graphState;                   % suppress graphs if necessary
    parameters=applyCorrectionWbreath(sourcePathName,fileList,destinationPathName,parameters,graphState);
    
    if parameters.Simulation.manualState    % reflect changed parameters in Calibration section
        fullfile(sourcePathName,cell2mat(fileList(end)));
        set(handles.settingsPath,   'String','settingPath');    % reset parameter settings file
    end
    
    set(handles.output,         'UserData',parameters);
    
catch err
    fclose('all');
    warning(sprintf('(batchProcessingButton_Callback): %s (in %s on line %g)\n', err.message,err.stack(1).name,err.stack(1).line))
end % try/catch

set(interfaceObject,'Enable','on');



% --- Executes on button press in graphButton.
function graphButton_Callback(hObject, eventdata, handles)
% hObject    handle to graphButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of graphButton
graphState=get(hObject,'Value');
parameters=get(handles.output,'UserData');
parameters.Simulation.graphState=graphState;
set(handles.output,'UserData',parameters);
fprintf('graphState is equal %d\n',graphState)



% --- Executes on button press in minVentButton.
function minVentButton_Callback(hObject, eventdata, handles)
% hObject    handle to minVentButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of minVentButton
minVentState=get(hObject,'Value');
parameters=get(handles.output,'UserData');
parameters.Simulation.minVentState=minVentState;
set(handles.output,'UserData',parameters);
fprintf('minute ventilation in logfile = %d\n',minVentState)



% --- Executes on button press in deleteLogButton.
function deleteLogButton_Callback(hObject, eventdata, handles)
% hObject    handle to deleteLogButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if exist('logFile.txt','file')
    answerString = questdlg('Are you sure to delete the log file?','Warning');
    if strcmp(answerString,'Yes')
        disp('deleting log file');
        delete 'logFile.txt'
        if exist('logFile.txt','file')
            warning('(deleteLogButton_Callback): cannot delete log file (probably open in another process)');
        end
    end
else
    fprintf('log file does not exist (aready deleted?)\n');
end



% --- Executes on button press in openLogButton.
function openLogButton_Callback(hObject, eventdata, handles)
% hObject    handle to openLogButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if exist('logFile.txt','file')
    data=fileread('logFile.txt');
    fprintf('copy log file into clipboard\n');
else
    data='log file does not exist';
    disp(data);
    data=strcat('<',data,'>');
end
clipboard('copy', data);



% --------------------------------------------------------------------
function interactiveMethod_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to interactiveMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
parameters=get(handles.output,'UserData');

interactiveYes=get(handles.interactiveYesButton,'Value');
parameters.Simulation.interactiveYes=interactiveYes;
fprintf('interactiveYes is equal %d\n',interactiveYes);

interactiveNo=get(handles.interactiveNoButton,'Value');
parameters.Simulation.interactiveNo=interactiveNo;
fprintf('interactiveNo is equal %d\n',interactiveNo);

interactiveDeviation=get(handles.interactiveDeviationButton,'Value');
parameters.Simulation.interactiveDeviation=interactiveDeviation;
fprintf('interactiveDeviation is equal %d\n',interactiveDeviation);

backgroundColor=get(0,'defaultUicontrolBackgroundColor');
foregroundColor=get(0,'defaultUicontrolForegroundColor');
if ~interactiveDeviation 
    set(handles.rmsCrit,'Background',backgroundColor,'Foreground',[0.5 0.5 0.5]);
else
    set(handles.rmsCrit,'Background','white','Foreground',foregroundColor);
end

set(handles.output,'UserData',parameters);



function rmsCrit_Callback(hObject, eventdata, handles)
% hObject    handle to rmsCrit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rmsCrit as text
%        str2double(get(hObject,'String')) returns contents of rmsCrit as a double
parameters=get(handles.output,'UserData');
parameters.Simulation.rmsCrit=str2double(get(hObject,'String'))/100;    % expressing as ratio instead of %
set(handles.output,'UserData',parameters);



% --- Executes during object creation, after setting all properties.
function rmsCrit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rmsCrit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in deleteAnonymousButton.
function deleteAnonymousButton_Callback(hObject, eventdata, handles)
% hObject    handle to deleteAnonymousButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

parameters=get(handles.output,'UserData');

[nEntries,mEntries]=size(parameters.Simulation.nameHash);
parameters.Simulation.nameHash = cell(0,2);
fprintf('hash table reset: %d entries deleted\n',nEntries);

set(handles.output,'UserData',parameters);



% --- Executes on button press in onlyPlotButton.
function onlyPlotButton_Callback(hObject, eventdata, handles)
% hObject    handle to onlyPlotButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of onlyPlotButton
parameters=get(handles.output,'UserData');
onlyPlotState=get(hObject,'Value');
parameters.Simulation.onlyPlotState=onlyPlotState;
fprintf('only plot state is equal %d\n',onlyPlotState)
set(handles.output,'UserData',parameters);



% --- Executes on button press in tidalButton.
function tidalButton_Callback(hObject, eventdata, handles)
% hObject    handle to tidalButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tidalButton
parameters=get(handles.output,'UserData');
parameters.Simulation.tidalState=get(hObject,'Value');
fprintf('Tidal state is equal %d\n',parameters.Simulation.tidalState)
if parameters.Simulation.tidalState || parameters.Simulation.croppingState
    set(handles.washoutButton,          'Enable', 'Off')
    set(handles.washinButton,           'Enable', 'Off')
else
    set(handles.washoutButton,          'Enable', 'On')
    set(handles.washinButton,           'Enable', 'On')    
end
if parameters.Simulation.tidalState
%     set(handles.minVentButton,          'Enable', 'Off')
    set(handles.lciStandardButton,      'Enable', 'Off')
    set(handles.lciByFitButton,         'Enable', 'Off')
    set(handles.TOevaluationButton,     'Enable', 'Off')
%    set(handles.momentRatiosButton,     'Enable', 'Off')
else
%     set(handles.minVentButton,          'Enable', 'On')
    set(handles.lciStandardButton,      'Enable', 'On')
    set(handles.lciByFitButton,         'Enable', 'On')
    set(handles.TOevaluationButton,     'Enable', 'On')
%    set(handles.momentRatiosButton,     'Enable', 'Off')
end
set(handles.output,'UserData',parameters);



% --- Executes on button press in croppingButton.
function croppingButton_Callback(hObject, eventdata, handles)
% hObject    handle to croppingButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of croppingButton
parameters=get(handles.output,'UserData');

croppingState=get(hObject,'Value');
parameters.Simulation.croppingState=croppingState;
fprintf('cropping state is equal %d\n',croppingState)
if parameters.Simulation.croppingState || parameters.Simulation.tidalState
    set(handles.washoutButton,          'Enable', 'Off')
    set(handles.washinButton,           'Enable', 'Off')
else
    set(handles.washoutButton,          'Enable', 'On')
    set(handles.washinButton,           'Enable', 'On')    
end
set(handles.output,'UserData',parameters);



% --- Executes on button press in lciStandardButton.
function lciStandardButton_Callback(hObject, eventdata, handles)
% hObject    handle to lciStandardButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lciStandardButton
parameters=get(handles.output,'UserData');

lciStandardState=get(hObject,'Value');
parameters.Simulation.lciStandard=lciStandardState;
fprintf('LCI standard = %d\n',lciStandardState);

set(handles.output,'UserData',parameters);



% --- Executes on button press in lciByFitButton.
function lciByFitButton_Callback(hObject, eventdata, handles)
% hObject    handle to lciByFitButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lciByFitButton
parameters=get(handles.output,'UserData');

lciByFitState=get(hObject,'Value');
parameters.Simulation.lciByFit=lciByFitState;
fprintf('LCI by fit = %d\n',lciByFitState)

set(handles.output,'UserData',parameters);



% --- Executes on button press in TOevaluationButton.
function TOevaluationButton_Callback(hObject, eventdata, handles)
% hObject    handle to TOevaluationButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TOevaluationButton
parameters=get(handles.output,'UserData');

TOevaluationState=get(hObject,'Value');
parameters.Simulation.TOevaluation=TOevaluationState;
fprintf('TO based evaluation = %d\n',TOevaluationState)

set(handles.output,'UserData',parameters);



% --- Executes on button press in momentRatiosButton.
function momentRatiosButton_Callback(hObject, eventdata, handles)
% hObject    handle to momentRatiosButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of momentRatiosButton
parameters=get(handles.output,'UserData');

momentRatiosState=get(hObject,'Value');
parameters.Simulation.momentRatios=momentRatiosState;
fprintf('moment ratios evaluation = %d\n',momentRatiosState);

set(handles.output,'UserData',parameters);



function vPre_Callback(hObject, eventdata, handles)
% hObject    handle to vPre (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vPre as text
%        str2double(get(hObject,'String')) returns contents of vPre as a double
parameters=get(handles.output,'UserData');
parameters.Device.volumePrecap=str2double(get(hObject,'String'));
set(handles.output,'UserData',parameters);



% --- Executes during object creation, after setting all properties.
function vPre_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vPre (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vPost_Callback(hObject, eventdata, handles)
% hObject    handle to vPost (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vPost as text
%        str2double(get(hObject,'String')) returns contents of vPost as a double
parameters=get(handles.output,'UserData');
parameters.Device.volumePostcap=str2double(get(hObject,'String'));
set(handles.output,'UserData',parameters);



% --- Executes during object creation, after setting all properties.
function vPost_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vPost (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in washoutButton.
function washoutButton_Callback(hObject, eventdata, handles)
% hObject    handle to washoutButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of washoutButton
parameters=get(handles.output,'UserData');
parameters.Simulation.washout=get(hObject,'Value');
fprintf('Washout state = %d\n',parameters.Simulation.washout)
set(handles.output,'UserData',parameters);



% --- Executes on button press in washinButton.
function washinButton_Callback(hObject, eventdata, handles)
% hObject    handle to washinButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of washinButton
parameters=get(handles.output,'UserData');
parameters.Simulation.washout=~get(hObject,'Value');
fprintf('Washin state = %d\n',~parameters.Simulation.washout)
set(handles.output,'UserData',parameters);



function lowerBoundEE_Callback(hObject, eventdata, handles)
% hObject    handle to lowerBoundEE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lowerBoundEE as text
%        str2double(get(hObject,'String')) returns contents of lowerBoundEE as a double
parameters=get(handles.output,'UserData');
parameters.Simulation.lowerBoundEE=str2double(get(hObject,'String'))/100;
set(handles.output,'UserData',parameters);



% --- Executes during object creation, after setting all properties.
function lowerBoundEE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lowerBoundEE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upperBoundEE_Callback(hObject, eventdata, handles)
% hObject    handle to upperBoundEE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of upperBoundEE as text
%        str2double(get(hObject,'String')) returns contents of upperBoundEE as a double
parameters=get(handles.output,'UserData');
parameters.Simulation.upperBoundEE=str2double(get(hObject,'String'))/100;
set(handles.output,'UserData',parameters);



% --- Executes during object creation, after setting all properties.
function upperBoundEE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to upperBoundEE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in meanButtonEE.
function meanButtonEE_Callback(hObject, eventdata, handles)
% hObject    handle to meanButtonEE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of meanButtonEE
parameters=get(handles.output,'UserData');
if get(hObject,'Value') && parameters.Simulation.meanEE~=1
else
    set(hObject,'Value',1)
end
parameters.Simulation.meanEE=1;
fprintf('MMee state = %d\n',parameters.Simulation.meanEE)
set(handles.output,'UserData',parameters);



% --- Executes on button press in medianButtonEE.
function medianButtonEE_Callback(hObject, eventdata, handles)
% hObject    handle to medianButtonEE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of medianButtonEE
parameters=get(handles.output,'UserData');
if get(hObject,'Value') && parameters.Simulation.meanEE~=2
else
    set(hObject,'Value',1)
end
parameters.Simulation.meanEE=2;
fprintf('MMee state = %d\n',parameters.Simulation.meanEE)
set(handles.output,'UserData',parameters);



% --- Executes on button press in byfitbuttonee.
function byFitButtonEE_Callback(hObject, eventdata, handles)
% hObject    handle to byfitbuttonee (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of byfitbuttonee
parameters=get(handles.output,'UserData');
parameters.Simulation.meanEE=~get(hObject,'Value');
if get(hObject,'Value') && parameters.Simulation.meanEE~=3
else
    set(hObject,'Value',1)
end
parameters.Simulation.meanEE=3;
fprintf('MMee state = %d\n',parameters.Simulation.meanEE)
set(handles.output,'UserData',parameters);



function nBreathMax_Callback(hObject, eventdata, handles)
% hObject    handle to nBreathMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nBreathMax as text
%        str2double(get(hObject,'String')) returns contents of nBreathMax as a double
parameters=get(handles.output,'UserData');
parameters.Simulation.nBreathMax=str2double(get(hObject,'String'));
set(handles.output,'UserData',parameters);



% --- Executes during object creation, after setting all properties.
function nBreathMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nBreathMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lowerBoundEI_Callback(hObject, eventdata, handles)
% hObject    handle to lowerBoundEI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lowerBoundEI as text
%        str2double(get(hObject,'String')) returns contents of lowerBoundEI as a double
parameters=get(handles.output,'UserData');
parameters.Simulation.lowerBoundEI=str2double(get(hObject,'String'))/100;
set(handles.output,'UserData',parameters);



% --- Executes during object creation, after setting all properties.
function lowerBoundEI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lowerBoundEI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upperBoundEI_Callback(hObject, eventdata, handles)
% hObject    handle to upperBoundEI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of upperBoundEI as text
%        str2double(get(hObject,'String')) returns contents of upperBoundEI as a double
parameters=get(handles.output,'UserData');
parameters.Simulation.upperBoundEI=str2double(get(hObject,'String'))/100;
set(handles.output,'UserData',parameters);



% --- Executes during object creation, after setting all properties.
function upperBoundEI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to upperBoundEI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in meanButtonEE.
function meanButtonEI_Callback(hObject, eventdata, handles)
% hObject    handle to meanButtonEE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of meanButtonEE
parameters=get(handles.output,'UserData');
if get(hObject,'Value') && parameters.Simulation.meanEI~=1
else
    set(hObject,'Value',1)
end
parameters.Simulation.meanEI=1;
fprintf('MMei state = %d\n',parameters.Simulation.meanEI)
set(handles.output,'UserData',parameters);



% --- Executes on button press in medianButtonEE.
function medianButtonEI_Callback(hObject, eventdata, handles)
% hObject    handle to medianButtonEE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of medianButtonEE
parameters=get(handles.output,'UserData');
if get(hObject,'Value') && parameters.Simulation.meanEI~=2
else
    set(hObject,'Value',1)
end
parameters.Simulation.meanEI=2;
fprintf('MMei state = %d\n',parameters.Simulation.meanEI)
set(handles.output,'UserData',parameters);



% --- Executes on button press in byfitbuttonee.
function byFitButtonEI_Callback(hObject, eventdata, handles)
% hObject    handle to byfitbuttonee (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of byfitbuttonee
parameters=get(handles.output,'UserData');
if get(hObject,'Value') && parameters.Simulation.meanEI~=3
else
    set(hObject,'Value',1)
end
parameters.Simulation.meanEI=3;
fprintf('MMei state = %d\n',parameters.Simulation.meanEI)
set(handles.output,'UserData',parameters);



% --- Executes on button press in fixedDeltaButton.
function fixedDeltaButton_Callback(hObject, eventdata, handles)
% hObject    handle to fixedDeltaButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fixedDeltaButton
parameters=get(handles.output,'UserData');
parameters.Simulation.deltaMMState=get(hObject,'Value');
fprintf('deltaMMFixed state is equal %d\n',parameters.Simulation.deltaMMState)
if parameters.Simulation.deltaMMState==1
	set(handles.deltaMMFix,'Enable','On')
else
	set(handles.deltaMMFix,'Enable','Off')
end
set(handles.output,'UserData',parameters);



function deltaMMFix_Callback(hObject, eventdata, handles)
% hObject    handle to deltaMMFix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of deltaMMFix as text
%        str2double(get(hObject,'String')) returns contents of deltaMMFix as a double
parameters=get(handles.output,'UserData');
parameters.Simulation.deltaMMFixed=str2double(get(hObject,'String'))/1000;
set(handles.output,'UserData',parameters);



% --- Executes during object creation, after setting all properties.
function deltaMMFix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to deltaMMFix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in tidalMeanButton.
function tidalMeanButton_Callback(hObject, eventdata, handles)
% hObject    handle to tidalMeanButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tidalMeanButton
parameters=get(handles.output,'UserData');

tidalMeanState=get(hObject,'Value');
parameters.Simulation.tidalMeanState=tidalMeanState;
fprintf('Tidal mean state = %d\n',tidalMeanState);

set(handles.output,'UserData',parameters);



function fromBreathTidal_Callback(hObject, eventdata, handles)
% hObject    handle to fromBreathTidal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fromBreathTidal as text
%        str2double(get(hObject,'String')) returns contents of fromBreathTidal as a double
parameters=get(handles.output,'UserData');
fromValue=str2double(get(hObject,'String'));
parameters.Simulation.fromBreathTidal=fromValue;
set(handles.output,'UserData',parameters);



% --- Executes during object creation, after setting all properties.
function fromBreathTidal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fromBreathTidal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nBreathTidal_Callback(hObject, eventdata, handles)
% hObject    handle to nBreathTidal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nBreathTidal as text
%        str2double(get(hObject,'String')) returns contents of nBreathTidal as a double
parameters=get(handles.output,'UserData');
nValue=str2double(get(hObject,'String'));
parameters.Simulation.nBreathTidal=nValue;
set(handles.output,'UserData',parameters);



% --- Executes during object creation, after setting all properties.
function nBreathTidal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nBreathTidal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in SF6Button.
function SF6Button_Callback(hObject, eventdata, handles)
% hObject    handle to SF6Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SF6Button
parameters=get(handles.output,'UserData');
parameters.Simulation.SF6=get(hObject,'Value');
fprintf('SF6 state is %d\n',parameters.Simulation.SF6)
set(handles.HeButton,'Value',~parameters.Simulation.SF6);
if parameters.Simulation.SF6 == 1
    set(handles.sideStreamButton,'Enable','off');
else
    set(handles.sideStreamButton,'Enable','on');
end
set(handles.output,'UserData',parameters);



% --- Executes on button press in HeButton.
function HeButton_Callback(hObject, eventdata, handles)
% hObject    handle to HeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of HeButton
parameters=get(handles.output,'UserData');
parameters.Simulation.SF6=~get(hObject,'Value');
fprintf('SF6 state is %d\n',parameters.Simulation.SF6)
set(handles.SF6Button,'Value',parameters.Simulation.SF6);
if parameters.Simulation.SF6 == 1
    set(handles.sideStreamButton,'Enable','off');
else
    set(handles.sideStreamButton,'Enable','on');
end
set(handles.output,'UserData',parameters);



% --- Executes on button press in sideStreamButton.
function sideStreamButton_Callback(hObject, eventdata, handles)
% hObject    handle to sideStreamButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sideStreamButton
parameters=get(handles.output,'UserData');
parameters.Simulation.sideStream=get(hObject,'Value');
fprintf('Side stream state is %d\n',parameters.Simulation.sideStream)
set(handles.output,'UserData',parameters);
