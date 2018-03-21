%   read_raw_adcp_bin.m

% program to read in binary file of raw ADCP data
clc
clear all
close all

% Set path to readin files
main_path = '/Users/marionsofiaalberty/MATLAB/Solomon_Sea/Moorings/Data/Raw/';

% addpath([matlabroot '/toolbox/matlab/RDADCP_oct06v0/']);
% addpath('/Applications/MATLAB74/toolbox/plots')

for moor_depl = 1:6
  %% Standard Cases (which require no special processing)
  switch moor_depl
    case 1
      % Solomon Strait M1 300 kHz
      % ----------------
      moor_name = 'Solomon_M1';
      ADCPdir = 'ADCP/';
      % WH 300 khz upward at 80 m
      Nrecords = 14428;
      DeployDepth = 80;
      SerNo = 5307;
      ADCPtype = 'RDI300kHz';
      RawADCPFile = 'M1-0080m_RDI-300kHz-5307_20140307.000';
      
    case 2
      % Solomon Strait M1 75 kHz
      % ----------------
      moor_name = 'Solomon_M1';
      ADCPdir = 'ADCP/';
      % WH75 khz down at 102 m
      Nrecords = 14430;
      DeployDepth = 102;
      SerNo = 3427;
      ADCPtype = 'RDI75kHz';
      RawADCPFile = 'M1-0102m_RDI-75kHz-3427_20140307.000';
      
    case 3
      % Solomon Strait M2b 75 kHz
      % ----------------
      moor_name = 'Solomon_M2b';
      ADCPdir = 'ADCP/';
      % WH75 khz upward at 400 m
      Nrecords = 14406;
      DeployDepth = 400;
      SerNo = 1066;
      ADCPtype = 'RDI75kHz';
      RawADCPFile = 'M2b-0400m_RDI_75kHz-1066_20140306.000';
      
    case 4
      % Solomon Strait M3 300 kHz
      % ----------------
      moor_name = 'Solomon_M3';
      ADCPdir = 'ADCP/';
      % WH 300 khz upward at 80 m
      Nrecords = 14382;
      DeployDepth = 80;
      SerNo = 12143;
      ADCPtype = 'RDI300kHz';
      RawADCPFile = 'M3-0080m_RDI-300kHz-12143_21040305.000';
      
    case 5
      % Vitiaz Middle 75kHz
      % --------------
      moor_name = 'VitiazMiddle';
      ADCPdir = 'ADCP/';
      % WH75 khz downward at 334 m
      Nrecords = 14407;
      DeployDepth = 334;
      SerNo = 16811;
      ADCPtype = 'RDI75kHz';
      RawADCPFile = 'VM-334m_RDI-75kHz-16811-Down_20140318.000';
      
    case 6
      % St Georges East 75kHz
      % --------------
      moor_name = 'StGeorgesEast';
      ADCPdir = 'ADCP/';
      % WH75 khz downward at 154 m
      Nrecords = 14445;
      DeployDepth = 154;
      SerNo = 16768;
      ADCPtype = 'RDI75kHz';
      RawADCPFile = 'SGE-0154m_RDI-75kHz-16768_20140313.000';
      
      %     case ???
      %         % Solomon Strait M2a 300 kHz
      %         % ----------------
      %         moor_name = 'Solomon_M2a';
      %         ADCPdir = 'ADCP/';
      %         % Flow Quest 300 kHz upward at 80 m
      %         Nrecords = 0;
      %         DeployDepth = 80;
      %         SerNo = nan;
      %         ADCPtype = 'FQ300kHz';
      %         RawADCPFile = '.000';
  end
  
  fprintf('Reading %s',RawADCPFile)
  
  % rdradcp(fname,num_ave_ensembles,num_total_ensembles);
  % if you don't specify the number of ensembles to average, defaults to 5
  % if you don't specify the total number of ensembles to read then it
  % tries to read past the end of the file.
  fname = ([main_path moor_name '/' ADCPdir RawADCPFile]);
  [adcp,config] = rdradcp(fname,1,Nrecords);
  
  OutFile = ([main_path moor_name '/' ADCPdir moor_name '_' ...
    num2str(DeployDepth) 'm_' ADCPtype '_' num2str(SerNo) '.mat']);
  
  % Save data
  save(OutFile,'adcp','config')
end

%% Non-Standard Cases (which require special processing)
for moor_depl = 7:10
  switch moor_depl
    case 7
      % Vitiaz Middle 75kHz
      % --------------
      % First bin is bad!
      % Nan out vel data
      moor_name = 'VitiazMiddle';
      ADCPdir = 'ADCP/';
      % WH75 khz upward at 332 m
      Nrecords = 27283; %tried to see if it goes past and stopped here
      DeployDepth = 332;
      SerNo = 8998;
      ADCPtype = 'RDI75kHz';
      RawADCPFile = 'VM-332m_RDI-75kHz-8998-Up_20140318.000';
      
      % Read in raw .000
      fprintf('Reading %s',RawADCPFile)
      fname = ([main_path moor_name '/' ADCPdir RawADCPFile]);
      [adcp,config] = rdradcp(fname,1,Nrecords);
      
      % nan out 1st bin in vel
      adcp.east_vel(1,:) = nan;
      adcp.north_vel(1,:) = nan;
      adcp.vert_vel(1,:) = nan;
      adcp.error_vel(1,:) = nan;
      
      OutFile = ([main_path moor_name '/' ADCPdir moor_name '_' ...
        num2str(DeployDepth) 'm_' ADCPtype '_' num2str(SerNo) '.mat']);
      
      % Save data
      save(OutFile,'adcp','config')
    case 8
      %% St Georges East 300 khZ
      % --------------
      % HAS A FALTY COMPASS!!!!
      % Use STGE down ADCP to fix
      moor_name = 'StGeorgesEast';
      ADCPdir = 'ADCP/';
      % WH300 khz upward at 152 m
      Nrecords = 28834;
      DeployDepth = 152;
      SerNo = 16832;
      ADCPtype = 'RDI300kHz';
      RawADCPFile = 'SGE-0152m_RDI-300kHz-16832_20140312.000';
      
      % Load StGE Down ADCP to use for compass correction
      % WH75 khz downward at 154 m
      DeployDepth_DOWN = 154;
      SerNo_DOWN = 16768;
      ADCPtype_DOWN = 'RDI75kHz';
      InFile = ([main_path moor_name '/' ADCPdir moor_name '_' ...
        num2str(DeployDepth_DOWN) 'm_' ADCPtype_DOWN '_' ...
        num2str(SerNo_DOWN) '.mat']);
      load(InFile,'adcp')
      adcp_DOWN = adcp;
      clear adcp
      
      % Read in raw .000
      fprintf('Reading %s',RawADCPFile)
      fname = ([main_path moor_name '/' ADCPdir RawADCPFile]);
      [adcp,config] = rdradcp(fname,1,Nrecords);
      
      % Correct the UP ADCP data
      adcp = fixUpCompass(adcp,adcp_DOWN,1);
      
      OutFile = ([main_path moor_name '/' ADCPdir moor_name '_' ...
        num2str(DeployDepth) 'm_' ADCPtype '_' num2str(SerNo) '.mat']);
      % Save data
      save(OutFile,'adcp','config')
      
    case 9
      %% St Georges West 75kHz and 300 kHz
      % --------------
      moor_name = 'StGeorgesWest';
      ADCPdir = 'ADCP/';
      % DOWN PRESSURE SENSOR IS BAD!!!!
      % Use STGW UP for pressue + offset
      % WH75 khz downward at 154 m
      Nrecords_DOWN = 14449;
      DeployDepth_DOWN = 154;
      SerNo_DOWN = 8866;
      ADCPtype_DOWN = 'RDI75kHz';
      RawADCPFile_DOWN = 'SGW-0154m_RDI-75kHz-8866_20140313.000';
      % UP HAS A FALTY COMPASS!!!!
      % Use STGW down ADCP to fix
      % WH300 khz upward at 152 m
      Nrecords_UP = 28834;
      DeployDepth_UP = 152;
      SerNo_UP = 16833;
      ADCPtype_UP = 'RDI300kHz';
      RawADCPFile_UP = 'SGW-152m_RDI-300kHz-16833_20140312.000';
      
      % Read in raw DOWN .000
      fprintf('Reading %s',RawADCPFile_DOWN)
      fname_DOWN = ([main_path moor_name '/' ADCPdir RawADCPFile_DOWN]);
      [adcp_DOWN,config_DOWN] = rdradcp(fname_DOWN,1,Nrecords_DOWN);
      
      % Read in raw UP .000
      fprintf('Reading %s',RawADCPFile_UP)
      fname_UP = ([main_path moor_name '/' ADCPdir RawADCPFile_UP]);
      [adcp_UP,config_UP] = rdradcp(fname_UP,1,Nrecords_UP);
      
      % Fix DOWN pressure
      % Define offset
      offset = 2.25; % instruments are separated by 2.25 m
      adcp_DOWN.depth = interp1(adcp_UP.mtime,adcp_UP.depth,...
        adcp_DOWN.mtime) + offset;
      adcp_DOWN.pressure = interp1(adcp_UP.mtime,adcp_UP.pressure,...
        adcp_DOWN.mtime) + (offset*1000);
      adcp_DOWN.pressure_std = interp1(adcp_UP.mtime,...
        adcp_UP.pressure_std,adcp_DOWN.mtime);
      
      % Save DOWN data
      adcp = adcp_DOWN;
      config = config_DOWN;
      OutFile_DOWN = ([main_path moor_name '/' ADCPdir moor_name '_' ...
        num2str(DeployDepth_DOWN) 'm_' ADCPtype_DOWN '_' ...
        num2str(SerNo_DOWN) '.mat']);
      save(OutFile_DOWN,'adcp','config')
      clear adcp config
      
      % Fix UP compass
      adcp_UP = fixUpCompass(adcp_UP,adcp_DOWN,2);
      
      % Save UP data
      adcp = adcp_UP;
      config = config_UP;
      OutFile_UP = ([main_path moor_name '/' ADCPdir moor_name '_' ...
        num2str(DeployDepth_UP) 'm_' ADCPtype_UP '_' num2str(SerNo_UP) ...
        '.mat']);
      save(OutFile_UP,'adcp','config')
      
    case 10
      %% Solomon Strait M3 75 kHz
      % ----------------
      % BROKEN INTO 2 FILES!!!
      % MAY ALSO FIX MAGNITUDE ISSUE HERE ...
      moor_name = 'Solomon_M3';
      ADCPdir = 'ADCP/';
      DeployDepth = 102;
      SerNo = 14215;
      ADCPtype = 'RDI75kHz';
      % WH 75 khz down at 102 m part 1
      Nrecords_1 = 6453;
      RawADCPFile_1 = 'M3-0102m_RDI-75kHz-14215_20140305-1.000';
      
      % WH 75 khz down at 102 m part 2
      Nrecords_2 = 4172;
      RawADCPFile_2 = 'M3-0102m_RDI-75kHz-14215_20140305-2.000';
      
      % Read in 1st raw .000
      fprintf('Reading %s',RawADCPFile_1)
      fname_1 = ([main_path moor_name '/' ADCPdir RawADCPFile_1]);
      [adcp_1,config] = rdradcp(fname_1,1,Nrecords_1);
      
      % Read in 2nd raw .000
      fprintf('Reading %s',RawADCPFile_2)
      fname_2 = ([main_path moor_name '/' ADCPdir RawADCPFile_2]);
      [adcp_2,~] = rdradcp(fname_2,1,Nrecords_2);
      
      % Concatinate data into one structure
      adcp = adcp_1;
      % Get all the fields in the data structure
      fields = fieldnames(adcp_1);
      for k = 3:numel(fields)
        % Get current field
        eval(sprintf('kField = adcp_1.%s;',fields{k}))
        % Check dimensions
        if ~ismatrix(kField)
          eval(sprintf('adcp.%s = cat(3,adcp_1.%s,adcp_2.%s);',...
            fields{k},fields{k},fields{k}))
        else
          eval(sprintf('adcp.%s = cat(2,adcp_1.%s,adcp_2.%s);',...
            fields{k},fields{k},fields{k}))
        end
      end
      
      % Save data
      OutFile = ([main_path moor_name '/' ADCPdir moor_name '_' ...
        num2str(DeployDepth) 'm_' ADCPtype '_' num2str(SerNo) '.mat']);
      save(OutFile,'adcp','config')
  end
end