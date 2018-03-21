% This function imports raw SBE56 files (*.asc) that have been uploaded from
% the SBE39 equiptment using the SeaBird Electronics SEATERM software
%
% INPUTS:
%   asc_filename: the filename of the SBE39 raw data file (*.asc)
%
% OUTPUTS:
%   times: matlab formatted date and time
%   temp: temperature records
%
% USAGE:
%   called by the gui function cleanSTR.m
%
% Author: Oliver Vetter November 2012
% adapted Feb 2013 to incorporate the 39 with pressure.

function [times, temp, pressure] = importSBE39(asc_filename)

fid = fopen(asc_filename);
tline = fgetl(fid);

% Ignore the header lines by running a quick string find for the last header
while ~feof(fid)
    % find the last line of the header flag for identifying the end of the head
    if isempty(strfind(tline, 'start sample number'));
    tline = fgetl(fid);
    else    
        break
    end
end
 
% for colobar
CData = 1:64; CData = (CData'*CData)/64;

% Question box
press = questdlg('Does this file contain pressure?', ...
                         'Pressure Question', ...
                         'Yes', 'No', 'No');
% Import boxes for 39 or 39 with pressure
if strcmp(press,'Yes');
    msgbox ('STR SBE39 with Pressure Data Importing, please wait','Seabird 39P','custom',CData, gray(64));
    [param, ~] = textscan(fid,'%f%f%s%s','Delimiter',',');
    % assign parameters
    temp = param{1}; pressure = param{2}; date = param{3}; time = param{4};
else strcmp(press,'No');
    msgbox ('STR SBE39 Temperature Data Importing, please wait','Seabird 39','custom',CData, hot(64));
    [param, ~] = textscan(fid,'%f%s%s','Delimiter',',');
    temp = param{1}; date = param{2}; time = param{3};
    pressure = 'No';
end

nn = datevec(time(1));
matTime = datenum(time,'HH:MM:SS') - datenum(nn(1),01,01);
matDate = datenum(date);
times = (matTime + matDate);

clear fid param