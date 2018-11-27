% mk_mooringTimeseries.m

% Interpolate and subsample all the data onto hourly and 20 m vertically
% spaced bins for each MoorSPICE mooring.

clear all; close all; clc
warning off

set(0,'defaultaxesfontsize',12,'defaultaxeslinewidth',0.7,...
  'defaultlinelinewidth',1,'defaultpatchlinewidth',0.7,...
  'defaultFigureColor','white')

%% Set script options
params.print = 0; % 1 = save figures
params.save = 0; % 1 = save data into 'mat' file

%% Set paths and directories
% Personnal paths:
% driveName='/Users/cyrilgermineaud/Documents/MATLAB/';
% addpath([driveName 'Routines_Cyril'])
driveName='/Users/marionsofiaalberty/MATLAB/Solomon_Sea/';

% Inputs/Outputs directories:
params.pathin =[driveName 'Moorings/Data/'];
params.pathout = [driveName 'Moorings/Data/Gridded/'];
params.pathfig = [driveName 'Moorings/Figures/Gridded/'];


%% Set parameters for all moorings
dz = 20;    % [m]
dt = 1/24;  % [days]

% Choose case
for i_moor = 1:9
  %% Set up each mooring case
  switch i_moor
    case 1
      % Solomon M1
      params.moor_name = 'Solomon_M1';
      params.channel = 'SolomonStrait';
      params.startgrid_mat = datenum(2012,07,22,01,00,00);
      params.startgrid = datestr(params.startgrid_mat);
      params.stopgrid_mat = datenum(2014,03,06,22,00,00);
      params.stopgrid = datestr(params.stopgrid_mat);
      params.gridDepth = 2050;
      params.sillDepth = 2525;
      params.bottomDepth = 2050;
      params.lat = -4-57.48/60;
      params.lon = 153+6.024/60;
      params.magdec = 6.64;
    case 2
      % Solomon M2a
      params.moor_name = 'Solomon_M2a';
      params.channel = 'SolomonStrait';
      params.startgrid_mat = datenum(2012,07,22,01,00,00);
      params.startgrid = datestr(params.startgrid_mat);
      params.stopgrid_mat = datenum(2014,03,05,07,00,00);
      params.stopgrid = datestr(params.stopgrid_mat);
      params.gridDepth = 2530;
      params.sillDepth = 2525;
      params.bottomDepth = 2559;
      params.lat = -5-9.853/60;
      params.lon = 153+16.864/60;
      params.magdec = 6.71;
    case 3
      % Solomon M2b
      params.moor_name = 'Solomon_M2b';
      params.channel = 'SolomonStrait';
      params.startgrid_mat = datenum(2012,07,23);
      params.startgrid = datestr(params.startgrid_mat);
      params.stopgrid_mat = datenum(2014,03,06,02,00,00);
      params.stopgrid = datestr(params.stopgrid_mat);
      params.gridDepth = 2530;
      params.sillDepth = 2525;
      params.bottomDepth = 2710;
      params.lat = -5-9.448/60;
      params.lon = 153+19.937/60;
      params.magdec = 6.71;
    case 4
      % Solomon M3
      params.moor_name = 'Solomon_M3';
      params.channel = 'SolomonStrait';
      params.startgrid_mat = datenum(2012,07,22,01,00,00);
      params.startgrid = datestr(params.startgrid_mat);
      params.stopgrid_mat = datenum(2014,03,04,20,00,00);
      params.stopgrid = datestr(params.stopgrid_mat);
      params.gridDepth = 2530;
      params.sillDepth = 2525;
      params.bottomDepth = 2627;
      params.lat = -5-8.507/60;
      params.lon = 154+17.975/60;
      params.magdec = 6.71;
    case 5
      % St. George's East
      params.moor_name = 'StGeorgesEast';
      params.channel = 'StGeorges';
      params.startgrid_mat = datenum(2012,07,20,08,00,00);
      params.startgrid = datestr(params.startgrid_mat);
      params.stopgrid_mat = datenum(2014,03,08);
      params.stopgrid = datestr(params.stopgrid_mat);
      params.gridDepth = 1433;
      params.sillDepth = 1400;
      params.bottomDepth = 1433;
      params.lat = -4-6.17/60;
      params.lon = 152+33.804/60;
      params.magdec = 6.39;
    case 6
      % St. George's West
      params.moor_name = 'StGeorgesWest';
      params.channel = 'StGeorges';
      params.startgrid_mat = datenum(2012,07,20,02,00,00);
      params.startgrid = datestr(params.startgrid_mat);
      params.stopgrid_mat = datenum(2014,03,07,20,00,00);
      params.stopgrid = datestr(params.stopgrid_mat);
      params.gridDepth = 1230;
      params.sillDepth = 1400;
      params.bottomDepth = 1243;
      params.lat = -4-6.82/60;
      params.lon = 152+31.116/60;
      params.breaktime = datenum(2013,01,04,15,10,01);
      params.fallDepth = 150;     % [m] Insts above this depth fell
      params.magdec = 6.39;
    case 7
      % Vitiaz East
      params.moor_name = 'VitiazEast';
      params.channel = 'Vitiaz';
      params.startgrid_mat = datenum(2012,07,28,12,00,00);
      params.startgrid = datestr(params.startgrid_mat);
      params.stopgrid_mat = datenum(2014,03,14,05,00,00);
      params.stopgrid = datestr(params.stopgrid_mat);
      params.gridDepth = 900;
      params.sillDepth = 1070;
      params.bottomDepth = 900;
      params.lat = -5-54.96/60;
      params.lon = 147+50.05/60;
      params.breaktime = datenum(2013,09,08,09,34,00);
      params.fallDepth = 600;     % [m] Insts above this depth fell
      params.magdec = 5.64;
    case 8
      % Vitiaz Middle
      params.moor_name = 'VitiazMiddle';
      params.channel = 'Vitiaz';
      params.startgrid_mat = datenum(2012,07,29,01,00,00);
      params.startgrid = datestr(params.startgrid_mat);
      params.stopgrid_mat = datenum(2014,03,14,02,00,00);
      params.stopgrid = datestr(params.stopgrid_mat);
      params.gridDepth = 1130;
      params.sillDepth = 1070;
      params.bottomDepth = 1130;
      params.lat = -5-56.65/60;
      params.lon = 147+46.68/60;
      params.breaktime = datenum(2012,10,22,15,31,00);
      params.fallDepth = 330;     % [m] Insts above this depth fell
      params.magdec = 5.64;
    case 9
      % Vitiaz West
      params.moor_name = 'VitiazWest';
      params.channel = 'Vitiaz';
      params.startgrid_mat = datenum(2012,07,28,07,00,00);
      params.startgrid = datestr(params.startgrid_mat);
      params.stopgrid_mat = datenum(2014,03,13,18,00,00);
      params.stopgrid = datestr(params.stopgrid_mat);
      params.gridDepth = 970;
      params.sillDepth = 1070;
      params.bottomDepth = 980;
      params.lat = -5-58.69/60;
      params.lon = 147+39.96/60;
      params.breaktime = datenum(2013,07,21,14,41,00);
      params.fallDepth = 600;     % [m] Insts above this depth fell
      params.magdec = 5.64;
  end
  
  %% Initialize intermediate data grid
  
  % Establish grid spacing for Time and Depth
  time = params.startgrid_mat:dt:params.stopgrid_mat;
  z = transpose(10:dz:params.gridDepth);
  
  % Initialize intermediate structures needed to make final T, S, U, and V
  % for T
  TT.pres = [];
  TT.pdep = [];
  TT.snum = {};
  TT.inst = {};
  TT.temp = [];
  % for S
  SS.pres = [];
  SS.pdep = [];
  SS.snum = {};
  SS.inst = {};
  SS.psal = [];
  % for U/V
  UV.pres = [];
  UV.pdep = [];
  UV.snum = {};
  UV.inst = {};
  UV.U = [];
  UV.V = [];
  
  
  %% Add available Level_1 data
  %  Read in files, filter higher frequency data, interpolate onto regular
  %  time, then add data and information to intermediate grids.
  
  % Make list of files in Level_1
  L1path = [params.pathin 'Level_1/' params.moor_name '/'];
  files = dir(fullfile(L1path,'*.mat'));
  totf = length(files);
  
  for k = 1:totf
    % Load file
    eval(['load ' L1path files(k).name])
    % Change structure name
    if exist('t','var')
      data=t;
    else
      data=s;
    end
    clear t s
    
    % Get instrument sampling frequency
    dt_inst = mean(diff(data.time));
    
    % If more frequently than twice hourly, filter
    if dt_inst < dt/2
      data = instFilter(data,dt);
    end
    
    % Interpolate data to regular time grid
    data = instTimeInterp(data,time);
    
    % Add data to intermediate grids
    [TT,SS,UV] = instAppend(data,TT,SS,UV);
  end
  
  
  %% Add remaining Level_0 data
  %  Read in file and serial number if not already in intermediate grid
  %  proceed with filtering higher frequency data, interpolation onto
  %  regular grid, and then add data and information to intermediate grids.
  
  % Make list of files in Level_0
  L0path = [params.pathin 'Level_0/' params.moor_name '/'];
  files = dir(fullfile(L0path,'*.mat'));
  totf = length(files);
  
  for k = 1:totf
    % Load file
    eval(['load ' L0path files(k).name])
    % Change structure name
    if exist('t','var')
      data=t;
    else
      data=s;
    end
    clear t s
    
    % Ignore instrument if a level_1 version was already used
    if ~SNexists(data,TT)
      
      % Get instrument sampling frequency
      dt_inst = mean(diff(data.time));
      
      % If more frequently than twice hourly, filter
      if dt_inst < dt/2
        data = instFilter(data,dt);
      end
      
      % Interpolate data to regular time grid
      data = instTimeInterp(data,time);
      
      % Add data to intermediate grids
      [TT,SS,UV] = instAppend(data,TT,SS,UV);
    end
  end
  clear data k
  
  % Ensure no 'above the surface' velocity data makes its way in
  UV.U(UV.pres < 0) = nan;
  UV.V(UV.pres < 0) = nan;
  
  
  %% Reorder data to be in order of planned meter depth
  % Sort UV
  UV = moorSort(UV);
  % Sort TT
  TT = moorSort(TT);
  % Sort SS
  SS = moorSort(SS);
  
  
  %% Deal with temporal offsets between instruments
  % Correct time offset using temp for TT
  TT = moorTimeOffset(TT,'temp',time);
  % Correct time offset using pres for SS
  SS = moorTimeOffset(SS,'pres',time);
  % Correct time offset using pres for UV
  UV = moorTimeOffset(UV,'pres',time);
  
  
  %% Interpolate pressure for instruments without pressure in temp struct
  % For moorings that had a break during deployment
  if ~isempty(TT.temp)
    if isfield(params,'breaktime')
      % Interpolate for pressure with break
      TT = moorFillPres_WB(TT,params,time);
    else
      % Interpolate for pressure for entire record
      TT = moorFillPres(TT);
    end
  end
  
  
  %% Apply the magnetic declination
  [UV.U,UV.V] = uvrot(UV.U,UV.V,params.magdec);
  
  
  %% Make plot of stacked observations
  if params.print
    moorPlotCoverage(TT,UV,time,params)
    disp('This is a chance to look at the data and stop if it looks bad.')
    disp('Press enter when you are done inspecting the plots.')
    pause
    disp('Here is your chance to stop if data looked bad.')
    prompt = 'Type "N" to stop program, anything else to continue.\n';
    user_choice = input(prompt,'s');
    if strcmp(user_choice,'N')
      return
    end
    close all
  end
  
  %% Interpolate (but not extrapolate) all data onto regular depth grid
  % Interpolate salinity to regular vertical grid
  if ~isempty(SS.psal)
    % Find T for RHO which only has the same sensors as S
    [~,~,iT] = intersect(SS.snum,TT.snum,'stable');
    % Make T for S for TS plots
    T4S.temp = TT.temp(iT,:);
    T4S.pres = TT.pres(iT,:);
    T4S.pdep = TT.pdep(iT);
    T4S.snum = TT.snum(iT);
    T4S.inst = TT.inst(iT);
    % Make RHO
    RHO.sgth = sw_pden(SS.psal,T4S.temp,SS.pres,0);
    RHO.pres = SS.pres;
    RHO.pdep = SS.pdep;
    RHO.snum = SS.snum;
    RHO.inst = SS.inst;
    % Interpolate
    SS = moorInterpZ(SS,z,time,'linear');
    RHO = moorInterpZ(RHO,z,time,'linear');
    T4S = moorInterpZ(T4S,z,time,'linear');
  end
  
  % Interpolate temperatue to regular vertical grid
  if ~isempty(TT.temp)
    TT = moorInterpZ(TT,z,time,'linear');
  end
  
  % Interpolate velocity to regular vertical grid
  UV = moorInterpZ(UV,z,time,'linear');
  
  
  %% Save data
  f_out = [params.pathout params.channel '/' params.moor_name '/' ...
    params.moor_name '.mat'];
  
  if ~isempty(SS.psal)
    save(f_out,'TT','SS','RHO','T4S','UV','params')
  elseif ~isempty(TT.temp)
    save(f_out,'TT','UV','params')
  else
    save(f_out,'UV','params')
  end
  
  % Clear variables
  clear UV TT SS RHO T4S
end