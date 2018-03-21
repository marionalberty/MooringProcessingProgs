% This function imports raw SBE56 files (*.csv) that have been uploaded from
% the SBE56 equiptment using the SeaBird Electronics SEATERM V2 software
%
% INPUTS:
%  -filename: the filename of the SBE56 raw data file (*.csv)
%  the csv file HAS to be exported without a header and ideally in a 
%  Windows Comma Separated foramat, in the
%   
%  SBE56 DATA NEED to be exported from SeaTermV2 'Without Headers'
%
% OUTPUTS:
%   times: matlab formatted date and time
%   temp: temperature records
%
% USAGE:
%   called by the gui function cleanSTR.m
%
% Author: Oliver Vetter November 2012

function [times, temp] = importSBE56(csv_filename)

fid = fopen(csv_filename);
tline = fgetl(fid);
 
%This if loop runs if there is only one line of header and breaks the if
%loop
if ~isempty(strfind(tline, '"Date","Time","Temperature"'));    
    [param, ~] = textscan(fid,'%s%s%f','Delimiter','","');
    % The "" quotation marks really mess up the dataread... this forces it
    % though
    data = param{2}; date = data(1:3:end); time = data(2:3:end);
    temp = data(3:3:end); temp = str2num(char(temp));
elseif ~isempty(strfind(tline, 'Date,Time,Temperature'));
    % Reads the file from line after the headerline 'start sample number = 1'
    [param, ~] = textscan(fid,'%s%s%f','Delimiter',',');
    % Assign parameters
    temp = param{3}; date = param{1}; time = param{2};
else errordlg ('!! Dataread Fail. Export SBE56 Without Headers !!','Export SBE56 Without Headers');
end

nn = datevec(time(1));
matTime = datenum(time,'HH:MM:SS') - datenum(nn(1),01,01);

% There are a number of ways to format the date from the software
matDate = datenum(date);
times = (matTime + matDate);

clear fid param