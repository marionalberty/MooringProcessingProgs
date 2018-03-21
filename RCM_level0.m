% RCM_level0

% Load raw.mat RCM files and calculate the associated speed and direction
% error for each time series. Save error in the structure as t.error.speed
% and t.error.direction.

clear all; close all; clc

set(0,'defaultaxesfontsize',12,'defaultaxeslinewidth',0.7,...
  'defaultlinelinewidth',1,'defaultpatchlinewidth',0.7,...
  'defaultFigureColor','white')


%% Set script options

params.print = 1; % 1 = save figures
params.save = 1; % 1 = save data into 'mat' file
params.debug = 0; % 1 = debug mode
params.qcLevel_in = 'Raw'; % define start qcLevel type
params.qcLevel_out = 'Level_0'; % define end qcLevel type


%% Set paths and directories

% Personnal paths:
% driveName='/Users/cyrilgermineaud/Documents/MATLAB/';
% addpath([driveName 'Routines_Cyril'])
driveName='/Users/marionsofiaalberty/MATLAB/Solomon_Sea/';

% Inputs/Outputs directories:
params.pathfile_in =[driveName 'Moorings/Data/' params.qcLevel_in '/'];
params.pathfile_out =[driveName 'Moorings/Data/' params.qcLevel_out '/'];
params.pathfig = [driveName 'Moorings/Figures/' params.qcLevel_out '/'];


%% Establish moorings to be analyzed

mooring_folder = {'Solomon_M1';'Solomon_M2b'; 'Solomon_M3'};
mooring_name = {'SolomonM1';'SolomonM2b'; 'SolomonM3'};
inst = 'RCM';


%% Load each file, calculate error and save Level 0 structure

% Go through by mooring
for i = 2:length(mooring_folder)
  pathin = [params.pathfile_in mooring_folder{i} '/' inst '/'];
  % Find raw.mat RCM files for this mooring
  files = dir(fullfile(pathin,'*.mat'));
  % Load files one at a time
  for k = 1:length(files)
    load([pathin files(k).name])
    % Make all variables single row vectors
    t.time = t.time(:)';
    t.pressure = t.pressure(:)';
    t.temperature = t.temperature(:)';
    t.u = t.u(:)';
    t.v = t.v(:)';
    t.speed = t.speed(:)';
    t.direction = t.direction(:)';
    
    % Add start/stop dates
    t.startdate_mat = t.time(1);
    t.stopdate_mat = t.time(end);
    t.startData = datestr(t.time(1));
    t.endData = datestr(t.time(end));
    
    %% Prompt user to see if there is drop-out
    % Plot
      plot(t.u)
      axis tight
      grid on;box on
      
      % prompt
      prompt = 'Is there velocity drop-out? y/n';
      DO_yn = input(prompt,'s');
      while strcmp(DO_yn,'y')
        % if yes
        zoom on;
        disp('Zoom into the area where drop-out occured. Hit any key when done.\n')
        pause();
        zoom off;
        disp('Select start and end point of a single drop out period.\n')
        [i_DO,~] = ginput(2);
        i_DO = round(i_DO);
        
        % change data to nan
        t.u(i_DO(1):i_DO(2)) = nan;
        t.v(i_DO(1):i_DO(2)) = nan;
        t.speed(i_DO(1):i_DO(2)) = nan;
        t.direction(i_DO(1):i_DO(2)) = nan;
        
        % prompt again
        prompt = 'Is there another velocity drop-out? y/n';
        DO_yn = input(prompt,'s');
      end
    
    
    %% Calculate speed error
    uv_err = ones(2,length(t.speed))*0.01;
    uv_err(2,:) = t.speed*0.02;
    t.error.speed = max(uv_err);
    
    % Calculate directional error
    dir_err = ones(size(t.direction))*5;
    dir_err(t.speed > 1) = 7.5;
    dir_err(t.speed < 0.05) = 7.5;
    t.error.direction = dir_err;
    
    % Add error readme
    t.error.readme = ['Error bounds are based on Aanderaa RCM7 manual '...
      'specifications for recording vector averaging. All error '...
      'estimates are +/-, with speed in m/s and direction in degrees.'];
    
    % Plot level_0 data
    if params.print
      figure('position',[0 0 800 800])
      set(gcf,'units','centimeters')
      set(gcf,'papersize',[8 8])
      set(gcf,'paperposition',[0,0,8,8])
      % Plot pressure
      subplot(411)
      plot(t.time,t.pressure,'k')
      axis tight ij
      datetick('x','mmm-yy','keeplimits')
      ylabel('Pressure [db]')
      title(sprintf('%s_%s_%s_%04i_0',mooring_name{i},t.meterType,...
        t.serialNo,t.plannedMeterDepth),'interpreter','None','fontsize',20)
      grid on;box on
      % Plot temperature
      subplot(412)
      plot(t.time,t.temperature,'r')
      axis tight
      datetick('x','mmm-yy','keeplimits')
      ylabel('Temperature [^\circC]')
      grid on;box on
      % Plot zonal velocity
      subplot(413)
      plot(t.time,t.u,'g')
      axis tight
      datetick('x','mmm-yy','keeplimits')
      ylabel('Zonal velocity,u [m/s]')
      grid on;box on
      % Plot zonal velocity
      subplot(414)
      plot(t.time,t.v,'g')
      axis tight
      datetick('x','mmm-yy','keeplimits')
      ylabel('Meridional velocity,v [m/s]')
      grid on;box on
      
      nfop=sprintf('%s_%s_%s_%04i_0',mooring_folder{i},t.meterType,...
        t.serialNo,t.plannedMeterDepth);
      % print png
      nfoppng=sprintf('%s%s/PNG/%s.png',params.pathfig,...
        mooring_folder{i},nfop);
      print(gcf,'-dpng','-r200',nfoppng);
      % print pdf
      nfoppdf=sprintf('%s%s/PDF/%s.png',params.pathfig,...
        mooring_folder{i},nfop);
      print(gcf,'-dpdf','-r200',nfoppdf);
      
      close all
    end
    
    % Save t
    fout = sprintf('%s%s/%s_%s_%s_%04i_0.mat',params.pathfile_out,...
      mooring_folder{i},mooring_name{i},t.meterType,t.serialNo,...
      t.plannedMeterDepth);
    save(fout,'t')
    
    % Clear t
    clear t
  end
end