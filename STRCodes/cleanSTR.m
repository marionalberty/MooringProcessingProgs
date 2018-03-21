function varargout = cleanSTR(varargin)
%% CLEANSTR.m MATLAB code for cleanSTR.fig
%
% FUNCTION:
%   imports raw STR files from the seabird SBE39 or SBE56 temperature
%   loggers.
%
% INPUTS:
%   Input the file following prompts: the filename of the SBE39 raw data file or (*.asc)
%
% OUTPUTS:
%   date and time as a vector format and temperature records in a comma
%   delimited text file in the standardized CRED Data Product (.cdp) format
%
% USAGE:
%   > cleanSTR
%
% CALLS:
%   cleanSTR.fig (gui)
%   selectDataPloy.m
%   importSBE39.m
%   importSBE56.m
%
%% Author: Oliver Vetter November 2012
% adapted Feb 19th 2013 to incorporate the 39 with pressure.

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cleanSTR_OpeningFcn, ...
                   'gui_OutputFcn',  @cleanSTR_OutputFcn, ...
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
% --- Executes just before cleanSTR is made visible.
function cleanSTR_OpeningFcn(hObject, ~, handles, varargin)

% Prompts to Load the data file
[filename, pathname] = uigetfile({'*.asc;*.csv'},'Please select an STR Data File to Import');
filenm = fullfile(pathname,filename);

% this is just for the colorbar 
CData = 1:64; CData = (CData'*CData)/64;
% loops that looks for a .asc file (assumes SBE39) or a .csv file (assumes
% a SBE56)
if strfind(filename,'.asc');
    [handles.times, handles.temp, handles.pressure] = importSBE39(filenm);
    STRType = 'SBE39';
else strfind(filename,'.csv');
    msgbox ('STR SBE56 data importing, please wait','Seabird 56','custom',CData,cool(64));
    [handles.times, handles.temp] = importSBE56(filenm);
    STRType = 'SBE56';
    handles.pressure = 'No';
end

if strcmp(handles.pressure,'No');
    press = 0;
else
    handles.pressureI = handles.pressure;
    press = 1;
end

figure(1);
if press == 1;
    subplot 211    
    plot(handles.times,handles.temp,'--k'); grid;
    % This hold is important for the select polygon to work
    hold on
    datetick('x','keeplimits');
    set(gca,'fontweight','bold');set(gca,'fontsize',15);
    title(['STR Type = ' STRType '     /     ' 'File = ' filename(1:end-13)]); 
    ylabel('Temperature (^oC)');
    subplot 212   
    plot(handles.times,handles.pressure,'--k'); grid;
    % This hold is important for the select polygon to work
    datetick('x','keeplimits');
    set(gca,'fontweight','bold');set(gca,'fontsize',15);
    xlabel('Date & Time'); ylabel('Pressure');
else
    plot(handles.times,handles.temp,'--k'); grid;
    % This hold is important for the select polygon to work
    hold on
    datetick('x','keeplimits');
    set(gca,'fontweight','bold');set(gca,'fontsize',15);
    title(['STR Type = ' STRType '     /     ' 'File = ' filename(1:end-13)]); 
    ylabel('Temperature (^oC)');
end

% Assign some data sturf
handles.timesI = handles.times; 
handles.tempI = handles.temp;
handles.pressureI = handles.pressure;
handles.press = press;
handles.filename = filename(1:end-13);
handles.filenm = filenm;
handles.STRType = STRType;

% Plots the data automatically
% Choose default command line output for cleanSTR
handles.output = hObject;
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = cleanSTR_OutputFcn(~, ~, handles)             

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in SelectData.
function SelectData_Callback(hObject, ~, handles)

figure(1);
clf
plot(handles.timesI,handles.tempI,'-+g'); grid;
hold on;datetick('x','keeplimits');
set(gca,'fontweight','bold'); set(gca,'fontsize',15);
title(['STR Type: ' handles.STRType '    /    ' 'File: ' (handles.filename)]); 
xlabel('Date & Time'); ylabel('Temperature (^oC)');

[I,~] = selectDataPoly(handles.timesI,handles.tempI);

qbox = questdlg ('Sure You Want to Remove These Data?','Remove Data?','Yes','No','Yes');
switch qbox
    case 'Yes'
    handles.timesI = handles.timesI(I); handles.tempI = handles.tempI(I); 
        if handles.press == 1;
        handles.pressureI = handles.pressureI(I);
        else
        end
    case 'No'
end
    
figure(1); hold off;
plot(handles.timesI,handles.tempI,'-+g'); grid;
hold on; datetick('x','keeplimits');
set(gca,'fontweight','bold'); set(gca,'fontsize',15);
title(['STR Type: ' handles.STRType '    /    ' 'File: ' (handles.filename)]); 
xlabel('Date & Time'); ylabel('Temperature (^oC)');
handles.output = hObject; guidata(hObject, handles);

% --- Executes on button press in RePlotData.
function RePlotData_Callback(hObject, ~, handles)
figure(1)
clf
hold off
if handles.press == 1;
    subplot 211    
    plot(handles.times,handles.temp,'--k'); grid;
    % This hold is important for the select polygon to work
    hold on
    datetick('x','keeplimits');
    set(gca,'fontweight','bold');set(gca,'fontsize',15);
    title(['STR Type = ' handles.STRType '     /     ' 'File = ' handles.filename]); 
    ylabel('Temperature (^oC)');
    subplot 212   
    plot(handles.times,handles.pressure,'--k'); grid;
    % This hold is important for the select polygon to work
    datetick('x','keeplimits');
    set(gca,'fontweight','bold');set(gca,'fontsize',15);
    xlabel('Date & Time'); ylabel('Pressure');
    handles.pressureI = handles.pressure;
else
    plot(handles.times,handles.temp,'--k'); grid;
    % This hold is important for the select polygon to work
    hold on
    datetick('x','keeplimits');
    set(gca,'fontweight','bold');set(gca,'fontsize',15);
    title(['STR Type = ' handles.STRType '     /     ' 'File = ' handles.filename]); 
    ylabel('Temperature (^oC)');
end% reassigns the original times and temps to reset the indicies

handles.timesI = handles.times;  handles.tempI = handles.temp;
handles.output = hObject; guidata(hObject, handles);

% --- Executes on button press in SaveFile.
function SaveFile_Callback(~, ~, handles)

SaveFileNm = handles.filenm(1:end-4);
[filename, pathname] = uiputfile({'*.cdp'},'Save data as .cdp CRED Data Product file',SaveFileNm);
file = ([pathname filename]);
msgbox ('Saving .cdp file,  please wait.','STR');

% parse out the dates and times as a vector to make the .cdp file
V = datevec(handles.timesI);
if handles.press == 1;
    cleanData = [V handles.tempI handles.pressureI];
else
    cleanData = [V handles.tempI];
end

% Save the text file to 7 significant digits
dlmwrite(file,cleanData,'delimiter','\t','precision',7);

% Closes the plot
close(figure(1)); 
% Closes the GUI window
close(figure(handles.figure1));
