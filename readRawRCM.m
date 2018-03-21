% readSolomonRCM

% Read raw RCM netcdf files in, generate and save mat file of the raw data
% to add to the raw folder for each Solomon Strait mooring


clear all; close all; clc

set(0,'defaultaxesfontsize',12,'defaultaxeslinewidth',0.7,...
  'defaultlinelinewidth',1,'defaultpatchlinewidth',0.7,...
  'defaultFigureColor','white')

%% Set script options

params.print = 1; % 1= save figures
params.save = 1; % 1= save raw data into 'mat' file
params.debug = 0; % 1= debug mode
params.qcLevel = 'Raw'; % define qcLevel type

%% Set paths and directories

% Personnal paths:
% driveName='/Users/cyrilgermineaud/Documents/MATLAB/';
% addpath([driveName 'Routines_Cyril'])
driveName='/Users/marionsofiaalberty/MATLAB/Solomon_Sea/';

% Inputs/Outputs directories:
params.pathfile =[driveName 'Moorings/Data/' params.qcLevel '/'];
params.pathfig = [driveName 'Moorings/Figures/' params.qcLevel '/'];

%% Establish moorings to be analyzed

mooring = {'Solomon_M1';'Solomon_M2b'; 'Solomon_M3'};
lat = [-4.9580; -5.1552; -5.1408];
lon = [153.1000; 153.3281; 154.2968];
wdp = [2050; 2710; 2617];
inst = 'RCM';
time_base=datenum(1901,1,15);

%% Load, format and save each file

% Go through by mooring
for i=1:length(mooring)
  pathin = [params.pathfile mooring{i} '/' inst '/'];
  % Find RCM files for this mooring
  files = dir(fullfile(pathin,'M*'));
  folds = files(logical([files(:).isdir]));
  files = files(~logical([files(:).isdir]));
  totf = length(files);
  
  % Read files one at a time
  for k=1:totf
    fname = [pathin files(k).name];
    finfo = strsplit(folds(k).name,'-');
    
    tnc_name=strrep(files(k).name(1:end-5),'_','');
    % Read in the 7 variables
    time = ncread(fname,[tnc_name 'T']);
    pres = ncread(fname,'PRESSION');
    temp = ncread(fname,'TEMP');
    sped = ncread(fname,'CURRENT_SPEED')./100;
    dirc = rad2deg(ncread(fname,'CURRENT_DIR'));
%     u = ncread(fname,'U')./100;
%     v = ncread(fname,'V')./100;
    [v,u]=pol2cart(deg2rad(dirc),sped);
    
    % Convert time to datenum
    time=time_base+time./24;    
    % Build structure
    t.latitude = lat(i);
    t.longitude = lon(i);
    t.plannedMeterDepth = str2double(tnc_name(3:end));
    t.waterDepth = wdp(i);
    t.meterType = inst;
    t.serialNo = finfo{3};
    t.units = {'time' 'Matlab date'; 'temperature' 'ITS-90, deg C';...
      'pressure' 'db'; 'velocity' 'm/s'; 'direction' 'degrees'};
    t.time = time;
    t.pressure = pres;
    t.temperature = temp;
    t.u = u;
    t.v = v;
    t.speed = sped;
    t.direction = dirc;
    t.startdate_mat = time(1);
    t.stopdate_mat = time(end);
    
    % Plot data and save figures
    % u & v
    subplot(211)
    plot(time,u,'b',[time(1) time(end)],[0 0],'r')
    datetick('x','mmm-yy','keeplimits')
    axis tight
    ylabel('U (m/s)')
    title(sprintf('%s %s %s Velocity (Mean Pressure) %8.4f',finfo{1},...
      finfo{3},inst,nanmean(pres)))
    subplot(212)
    plot(time,v,'b',[time(1) time(end)],[0 0],'r')
    datetick('x','mmm-yy','keeplimits')
    axis tight
    ylabel('V (m/s)')
    nfoppng = sprintf('%s%s/%s_%sm_%s_RCM_Velocity.png',params.pathfig,...
      mooring{i},mooring{i},finfo{2},finfo{3});
    print(gcf,'-dpng',nfoppng);
    close
    
    % speed
    subplot(211)
    plot(time,sped,'b')
    datetick('x','mmm-yy','keeplimits')
    axis tight
    ylabel('Speed (m/s)')
    title(sprintf('%s %s %s Velocity (Mean Pressure) %8.4f',finfo{1},...
      finfo{3},inst,nanmean(pres)))
    subplot(212)
    plot(time,dirc,'b')
    datetick('x','mmm-yy','keeplimits')
    axis tight
    ylabel('Direction (degrees)')
    nfoppng = sprintf('%s%s/%s_%sm_%s_RCM_speed.png',params.pathfig,...
      mooring{i},mooring{i},finfo{2},finfo{3});
    print(gcf,'-dpng',nfoppng);
    close
    
    % t & p
    subplot(211)
    plot(time,pres,'b')
    datetick('x','mmm-yy','keeplimits')
    axis tight
    ylabel('Pressure (db)')
    title(sprintf('%s %s %s Velocity (Mean Pressure) %8.4f',finfo{1},...
      finfo{3},inst,nanmean(pres)))
    subplot(212)
    plot(time,temp,'b')
    datetick('x','mmm-yy','keeplimits')
    axis tight
    ylabel('Temperature (deg C)')
    nfoppng = sprintf('%s%s/%s_%sm_%s_RCM_Velocity_Pressure_T.png',...
      params.pathfig,mooring{i},mooring{i},finfo{2},finfo{3});
    print(gcf,'-dpng',nfoppng);
    close
    
    % Save t
    fout = sprintf('%s%s/RCM/%s_%s_%s_RCM_raw.mat',params.pathfile,...
      mooring{i},mooring{i},finfo{2},finfo{3});
    save(fout,'t')
    
    % Clear t
    clear t
  end
end


