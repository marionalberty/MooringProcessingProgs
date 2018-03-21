% adcp_qc.m

% Programmer: Janet Sprintall
% Date: 30 November 2007
% Edited and adapted by Marion Alberty for use on Solomon Sea moorings in
% June 2015
% Modified by Cyril Germineaud, put vars and params into structures and
% define plot_adcp_qc.m function
% Date: August 2015

% Purpose: program reads in the mat file of the raw adcp data
% and quality control the data according to NDBC QC document modification
% for the WH300 and the LR75 instruments.
% plots various diagnostic parameters from the QC procedure and the mooring
% (e.g. EA, pitch roll etc).

% input:
% Mcases        switch for each ADCP deployment to qc

% output to:
% adcp_file =   [main_path filein moor_name SerNo '_0.mat'];
% main_path     main directory
% filein        subdirectory for that mooring
% moor_name     mooring name and deployment number
%               Dipolog_deploy1
% SerNo         serial number

close all; clear all; clc

set(0,'defaultaxesfontsize',14,'defaultaxeslinewidth',0.7,...
  'defaultlinelinewidth',1,'defaultpatchlinewidth',0.7,...
  'defaultFigureColor','white')

% Personnal paths:
driveName='/Users/marionsofiaalberty/MATLAB/Solomon_Sea/';
% driveName='/Users/cyrilgermineaud/Documents/MATLAB/';
% addpath([driveName 'Routines_Cyril'])
% addpath([driveName 'ToolBox'])
% addpath([driveName 'ToolBox/graphics'])
% addpath([driveName 'ToolBox/graphics/paruly'])

% Set up paths
params.main_path = [driveName 'Moorings/'];
% params.main_path = [driveName 'Solomon_Sea_Moorings/'];
params.file_dat = 'Data/';

% path for figures
params.figs_raw = [params.main_path 'Figures/Raw/QC_tests/'];
params.figs_lev0 = [params.main_path 'Figures/Level_0/'];

% Set up script options:
params.print = 1;               % 1= print figures
params.save = 1;                % 1= save level0 data into 'mat' file
params.qcLevel = 0;             % define qcLevel type

% Mooring cases to QC
Mcases = 1:11;

%% 1. INITIALIZE INPUT FOR EACH MOORING

for idmoor = Mcases
  disp 'on going mooring instrument...'
  disp(idmoor)
  
  switch idmoor
    case 1
      % St Georges East 300 khZ
      % --------------
      % file name
      params.moor_name = 'StGeorgesEast';
      
      % WH300 khz upward at 152 m
      WaterDepth = 1433;
      DeployDepth = 152;
      SerNo = 16832;
      ADCPtype = 'RDI300';
      
      % adcp orientation
      aorient = 'up';
      lat = -4-6.174/60;
      lon = 152+31.116/60;
      
      % deployment times
      %start_time=datenum(2012,07,20,06,10,00);
      %end_time = datenum(2014,03,8,04,00,00);
      start_time=datenum(2012,07,20,07,04,24); % start data
      end_time = datenum(2014,03,8,04,34,24); % end data
      clock_drift = -11/3600; % convert into hrs
      
      % Magnetic declination not added to raw ADCP data:
      % add to structure and in processing
      % At 4deg6.2S, 152deg33.8E the magnetic declination is 6.39degE ...
      % (i.e. positive, move it clockwise)
      magdec = 6.39;
            
    case 2
      % St Georges East 75kHz
      % --------------
      % file name
      params.moor_name = 'StGeorgesEast';
      % LR75 khz downward at 152 m
      
      WaterDepth = 1433;
      DeployDepth = 154;
      SerNo = 16768;
      ADCPtype = 'RDI75';
      
      % adcp orientation
      aorient = 'down';
      lat = -4-6.174/60;
      lon = 152+31.116/60;
      
      % deployment times
      %start_time=datenum(2012,07,20,06,10,00);
      %end_time = datenum(2014,03,8,04,00,00);
      start_time=datenum(2012,07,20,06,20,00); %start data
      end_time = datenum(2014,03,8,03,20,00); %end data
      clock_drift = 100/3600; % convert into hrs
      
      % Magnetic declination
      magdec = 6.39;
      
    case 3
      % St Georges West 75kHz
      % --------------
      % file name
      params.moor_name = 'StGeorgesWest';
      % LR75 khz downward at 152 m
      
      WaterDepth = 1243;
      DeployDepth = 154;
      SerNo = 8866;
      ADCPtype = 'RDI75';
      
      % adcp orientation
      aorient = 'down';
      lat = -4-6.174/60;
      lon = 152+33.804/60;
      
      % deployment times
      %start_time=datenum(2012,07,20,06,00,00);
      %end_time = datenum(2014,03,7,21,00,00);
      start_time=datenum(2012,07,20,02,11,12); % start data
      end_time = datenum(2014,03,7,21,11,11);% end data
      clock_drift = -28/3600; % convert into hrs
      % Magnetic declination
      magdec = 6.39;
      
    case 4
      % St Georges West 300kHz
      % --------------
      % file name
      params.moor_name = 'StGeorgesWest';
      % LR75 khz downward at 152 m
      
      WaterDepth = 1243;
      DeployDepth = 152;
      SerNo = 16833;
      ADCPtype = 'RDI300';
      
      % adcp orientation
      aorient = 'up';
      lat = -4-6.174/60;
      lon = 152+33.804/60;
      
      % deployment times
      %start_time=datenum(2012,07,20,06,00,00);
      %end_time = datenum(2014,03,7,21,00,00);
      start_time=datenum(2012,07,20,01,32,24); % start data
      end_time = datenum(2014,03,7,20,32,24); % end data
      clock_drift = 69/3600; % convert into hrs
      
      % Magnetic declination
      magdec = 6.39;
      
    case 5
      % Vitiaz Middle 75kHz UP
      % --------------
      % file name
      params.moor_name = 'VitiazMiddle';
      
      % WH75 khz upward at 332 m
      WaterDepth = 1131;
      DeployDepth = 332;
      SerNo = 8998;
      ADCPtype = 'RDI75';
      
      % adcp orientation
      aorient = 'up';
      lat = -5-56.65/60;
      lon = 147+46.68/60;
      
      % deployment times
      %start_time=datenum(2012,07,28,06,15,00);
      %end_time = datenum(2014,03,13,18,45,00);
      start_time=datenum(2012,07,28,23,00,00); % start data
      end_time = datenum(2014,02,14,10,59,59); % end data
      clock_drift = 375/3600; % convert into hrs
      
      % Magnetic declination
      magdec = 5.64;
      
    case 6
      % Vitiaz Middle 75kHz DOWN
      % --------------
      % file name
      params.moor_name = 'VitiazMiddle';
      
      % WH75 khz downward at 334 m
      WaterDepth = 1131;
      DeployDepth = 334;
      SerNo = 16811;
      ADCPtype = 'RDI75';
      
      % adcp orientation
      aorient = 'down';
      lat = -5-56.65/60;
      lon = 147+46.68/60;
      
      % deployment times
      %start_time=datenum(2012,07,28,06,15,00);
      %end_time = datenum(2014,03,13,18,45,00);
      start_time=datenum(2012,07,29,00,10,24); % start data
      end_time = datenum(2014,03,14,01,10,24); % end data
      clock_drift = 124/3600; % convert into hrs
      
      % Magnetic declination
      magdec = 5.64;
      
    case 7
      % Solomon Strait M1 300 kHz UP
      % ----------------
      % file name
      params.moor_name = 'Solomon_M1';
      
      WaterDepth = 2100;
      DeployDepth = 80;
      SerNo = 5307;
      ADCPtype = 'RDI300';
      
      % adcp orientation
      aorient = 'up';
      lat = -4-57.48/60;
      lon = 153+6.024/60;
      
      % deployment times
      %start_time = datenum(2012,07,21,01,30,00);
      %end_time = datenum(2014,03,06,23,01,00);
      start_time = datenum(2012,07,20,23,24,00); % start data
      end_time = datenum(2014,03,06,19,24,00); % end data
      clock_drift = 540/3600; % convert into hrs
      
      % Magnetic declination
      magdec = 6.64;
      
    case 8
      % Solomon Strait M1 75 kHz DOWN
      % ----------------
      % file name
      params.moor_name = 'Solomon_M1';
      
      WaterDepth = 2100;
      DeployDepth = 102;
      SerNo = 3427;
      ADCPtype = 'RDI75';
      
      % adcp orientation
      aorient = 'down';
      lat = -4-57.48/60;
      lon = 153+6.024/60;
      
      % deployment times
      %start_time = datenum(2012,07,21,01,30,00);
      %end_time = datenum(2014,03,06,23,01,00);
      start_time = datenum(2012,07,20,21,28,00); % start data
      end_time = datenum(2014,03,06,17,28,00); % end data
      clock_drift = 830/3600; % convert into hrs
      
      % Magnetic declination
      magdec = 6.64;
      
    case 9
      % Solomon Strait M2b 75 kHz UP
      % ----------------
      % file name
      params.moor_name = 'Solomon_M2b';
      
      WaterDepth = 2710;
      DeployDepth = 400;
      SerNo = 1066;
      ADCPtype = 'RDI75';
      
      % adcp orientation
      aorient = 'up';
      lat = -5-9.448/60;
      lon = 153+19.937/60;
      
      % deployment times
      %start_time = datenum(2012,07,15,22,11,00);
      %end_time = datenum(2014,03,06,02,24,00);
      start_time = datenum(2012,07,15,21,37,36); % start data
      end_time = datenum(2014,03,06,00,37,36); % end data
      clock_drift = 206/3600; % convert into hrs
      
      % Magnetic declination
      magdec = 6.71;
      
    case 10
      % Solomon Strait M3 300 kHz UP
      % ----------------
      % file name
      params.moor_name = 'Solomon_M3';
      
      WaterDepth = 2617;
      DeployDepth = 80;
      SerNo = 12143;
      ADCPtype = 'RDI300';
      
      % adcp orientation
      aorient = 'up';
      lat = -5-8.507/60;
      lon = 154+17.975/60;
      
      % deployment times
      %start_time = datenum(2012,07,15,22,11,00);
      %end_time = datenum(2014,03,06,02,24,00);
      start_time = datenum(2012,07,17,00,50,24); % start data
      end_time = datenum(2014,03,04,15,50,24); % end data
      clock_drift = 774/3600; % convert into hrs
      
      % Magnetic declination
      magdec = 6.71;
      
    case 11
      % Solomon Strait M3 75 kHz DOWN
      % -----------------------------
      % file name
      params.moor_name = 'Solomon_M3';
      
      WaterDepth = 2617;
      DeployDepth = 102;
      SerNo = 14215;
      ADCPtype = 'RDI75';
      
      % adcp orientation
      aorient = 'down';
      lat = -5-9.448/60;
      lon = 153+19.937/60;
      
      % deployment times
      %start_time = datenum(2012,07,15,22,11,00);
      %end_time = datenum(2014,03,06,02,24,00);
      start_time = datenum(2012,07,16,23,14,48); % start data
      end_time = datenum(2013,09,29,23,18,11); % end data
      clock_drift = 563/3600; % convert into hrs
      
      % Magnetic declination
      magdec = 6.71;
      
      
      %     case ???
      %
      %       % Solomon Strait M2a 300 kHz UP
      %       % ----------------
      %       % for file names
      %       params.moor_name = 'Solomon_M2a';
      %
      %       WaterDepth = 2900;
      %       DeployDepth = 80;
      %       SerNo = 40005;
      %       ADCPtype = 'FQ300';
      %
      %       % adcp orientation
      %       aorient = 'up';
      %       lat = -5-9.853/60;
      %       lon = 153+16.864/60;
      %
      %       % deployment times
      %       start_time = datenum(2012,07,16,04,17,00);
      %       end_time = datenum(2014,03,06,01,45,00);
      %
      %       % Magnetic declination
      %       magdec = 6.70;
  end
  ADCPdir = 'ADCP/';
  
  
  %% 2. READ IN RAW DATA FILE AND SET UP THRESHOLD FOR ADCP QC TESTS
  
  filein = ([params.main_path params.file_dat 'Raw/' params.moor_name ...
    '/' ADCPdir params.moor_name '_' num2str(DeployDepth) 'm_' ADCPtype ...
    'kHz_' num2str(SerNo) '.mat']);
  
  % Load filein
  load(filein)
  
  % Get QC threshholds
  qcthresh = params_threshold_qc(params,adcp.config.beam_freq,...
    adcp.config.orientation);
  
  fprintf('Threshold Parameters for ADCP quality control\n\n')
  fprintf('Error Velocity %5.2f meters per second\n', qcthresh.err_vel)
  fprintf('Percent Good %5.2f \n',qcthresh.pgood)
  fprintf('Correlation Magnitude %5.2f \n',qcthresh.cmag)
  fprintf('Vertical Velocity %5.2f meters per second\n',qcthresh.vvel)
  fprintf('Speed %5.2f meters per second\n',qcthresh.hvel)
  fprintf('Echo Amplitude bin difference %5.2f \n\n',qcthresh.ea_thresh)
  
  % set up a few parameters that are needed
  time0=adcp.mtime;
  adcp.mtime = time0-clock_drift;
  nbins = adcp.config.n_cells;
  samplerate = (adcp.mtime(10)-adcp.mtime(9))*60*24; % in minutes
  uu = adcp.east_vel;
  vv = adcp.north_vel;
  u = uu + 1i* vv;
  w = adcp.vert_vel;
  erv = adcp.error_vel;
  
  % set up a structure qc for diagnostics
  for iz = 1:nbins
    for nn =1:4
      qc(nn).ea(iz,:) = adcp.intens(iz,nn,:);
      qc(nn).pg(iz,:) = adcp.perc_good(iz,nn,:);
      qc(nn).cr(iz,:) = adcp.corr(iz,nn,:);
    end
  end
  
  
  %% 3. RUN DIAGNOSTICS FOR ADCP QC CHECK
  
  disp('calling adcp qc test')
  
  % filen is the filename
  %[ifailar,ifailac] where >2 failed the first 5 tests
  %[ifailr,ifailc] failed the ea test
  
  % Define path for RAW figures per mooring
  params.subdir_raw = [params.figs_raw params.moor_name '/'];
  
  [ifailar,ifailac,ifailr,ifailc] = adcpqctest_2(qcthresh,qc',u,w,erv,...
    params.moor_name,ADCPtype,SerNo,params.subdir_raw);
  
  % set the velocity data to nan if failed qc test
  % first 5 tests
  
  for ir = 1:length(ifailar)
    u(ifailar(ir),ifailac(ir)) = complex(nan,nan);
    w(ifailar(ir),ifailac(ir)) = nan;
    erv(ifailar(ir),ifailac(ir)) = nan;
  end
  % ea tests
  for ir = 1:length(ifailr)
    u(ifailr(ir):end,ifailc(ir)) = complex(nan,nan);
    w(ifailr(ir):end,ifailc(ir)) = nan;
    erv(ifailr(ir):end,ifailc(ir)) = nan;
  end
  
  %% Now put all data onto common time base dept of start and end time
  
  isub = find(adcp.mtime >= start_time & adcp.mtime <= end_time);
  % original velocities
  uu = uu(:,isub);
  vv = vv(:,isub);
  v = imag(u(:,isub));
  u = real(u(:,isub));
  w = w(:,isub);
  erv = erv(:,isub);
  temperature = adcp.temperature(isub);
  pressure = adcp.pressure(isub)/1000;
  ndate = adcp.mtime(isub);
  depth = adcp.depth(isub);
  pitch = adcp.pitch(isub);
  pitch_std = adcp.pitch_std(isub);
  roll = adcp.roll(isub);
  roll_std = adcp.roll_std(isub);
  heading = adcp.heading(isub);
  heading_std = adcp.heading_std(isub);
  
  for ii =1:4
    qc(ii).ea = qc(ii).ea(:,isub);
    qc(ii).pg = qc(ii).pg(:,isub);
    qc(ii).cr = qc(ii).cr(:,isub);
  end
  
  % set up depth structure
  % test for depth record
  iz = find(adcp.depth == 0);
  if(length(iz)==length(adcp.depth))
    adcp.depth = nan(1,length(adcp.depth));
    adcp.pressure = adcp.depth;
    disp('no depth record for this ADCP')
  end
  
  % test for pressure record
  % if none then work out
  izp = find(adcp.pressure == 0);
  if(length(izp)==length(adcp.pressure))
    pressure = sw_pres(adcp.depth',lat);
    disp('no pressure record for this ADCP')
  end
  
  if strcmp(aorient(1:2),'up')
    depth = adcp.depth' * ones(1,length(adcp.config.ranges)) - ...
      ones(length(adcp.depth),1)*adcp.config.ranges';
  elseif strcmp(aorient(1:2),'do')
    depth = adcp.depth' * ones(1,length(adcp.config.ranges)) + ...
      ones(length(adcp.depth),1)*adcp.config.ranges';
  end
  
  depth = depth(isub,:);
  
  % reject negative depth bins
  gibad = find(depth' < 0);
  u(gibad) = NaN;
  v(gibad) = NaN;
  w(gibad) = NaN;
  erv(gibad) = NaN;
  depth(depth < 0) = NaN;
  
  %% save data in structure
  % all q.c. adcp data
  s = struct('u',u);
  s = setfield(s,'v',v);
  s = setfield(s,'w',w);
  s = setfield(s,'erv',erv);
  s = setfield(s,'temperature',temperature);
  s = setfield(s,'pressure',pressure);
  s = setfield(s,'time',ndate);
  s = setfield(s,'depth',depth');
  s = setfield(s,'qcthresh',qcthresh);
  s = setfield(s,'magneticDeclination',magdec);
  s = setfield(s,'latitude',lat);
  s = setfield(s,'longitude',lon);
  s = setfield(s,'meterType',ADCPtype);
  s = setfield(s,'serialNo',num2str(SerNo));
  s = setfield(s,'waterDepth',WaterDepth);
  s = setfield(s,'plannedMeterDepth',DeployDepth);
  s = setfield(s,'sampleInterval',samplerate);
  s = setfield(s,'range',adcp.config.ranges);
  s = setfield(s,'mooringName',params.moor_name);
  s = setfield(s,'QCvariables',qc);
  s = setfield(s,'orientation',aorient);
  
  % add ndate, heading, pitch and roll into 's' struct
  s.startData = datestr(ndate(1));
  s.endData = datestr(ndate(end));
  s.heading = heading;
  s.pitch = pitch;
  s.roll = roll;
  s.heading_std = heading_std;
  s.pitch_std = pitch_std;
  s.roll_std = roll_std;
  
  params.file_out = sprintf('%s_%s_%s_%4.4d_%d.mat',...
    strrep(params.moor_name,'_',''),s.meterType,s.serialNo,...
    s.plannedMeterDepth,params.qcLevel);
  
  % Save ADCP data into the Level_0 folder
  if params.save==1
    eval(['save ' params.main_path params.file_dat 'Level_0/' ...
      params.moor_name '/' params.file_out ' s'])
  end
  
  %% 6. PLOT THE ADCP DATA
  
  close all
  
  % Define figure subdirectories per mooring
  params.dirpng = [params.figs_lev0 params.moor_name '/PNG/ADCP_QC'];
  params.dirpdf = [params.figs_lev0 params.moor_name '/PDF/ADCP_QC'];
  
  % Plot level 0 data and qc diagnostics
  plot_adcp_qc(params,nbins,s,uu,vv,qc)
  
  % Plot level 0 data
  plot_adcp_level0(params,s)
  
  % Clear variables
  clearvars -except driveName params
end