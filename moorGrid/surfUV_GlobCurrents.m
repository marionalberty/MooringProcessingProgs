% surfUV_GlobCurrents.m

% Load, interpolate, and save surface velocities using GlobCurrents for the
% Solomon Strait moorings.

clear all; close all; clc
warning off

set(0,'defaultaxesfontsize',16,'defaultaxeslinewidth',0.7,...
  'defaultlinelinewidth',1,'defaultpatchlinewidth',0.7,...
  'defaultFigureColor','white')


%% Set paths and directories
% Personnal paths:
% driveName='/Users/cyrilgermineaud/Documents/MATLAB/';
% addpath([driveName 'Routines_Cyril'])
driveName = '/Users/marionsofiaalberty/MATLAB/';
dataPath = 'Solomon_Sea/Moorings/Data/Gridded/';
gcurPath = 'GlobCurrent/MoorSPICE/';


% Choose cases
for i_moor = 1:3
  %% Set up each mooring case and load data
  switch i_moor
    case 1
      % Solomon M1
      params.moor_name = 'Solomon_M1';
      params.channel = 'SolomonStrait';
      load([driveName dataPath params.channel '/' params.moor_name ...
        '/Solomon_M1.mat'],'UV','params')
      params.lon = params.lon + 0.1;
    case 2
      % Solomon M2b
      params.moor_name = 'Solomon_M2b';
      params.channel = 'SolomonStrait';
      load([driveName dataPath params.channel '/' params.moor_name ...
        '/Solomon_M2b.mat'],'UV','params')
    case 3
      % Solomon M3
      params.moor_name = 'Solomon_M3';
      params.channel = 'SolomonStrait';
      load([driveName dataPath params.channel '/' params.moor_name ...
        '/Solomon_M3.mat'],'UV','params')
  end
  
  
  %% Load GlobCurrent surface velocities, interpolate
  load([driveName gcurPath ...
    'GLOBCURRENT_total_MoorSPICE_20120101-20141231_L4_v3.mat'],...
    'lat','lon','time','u_surface','v_surface')
  load([driveName gcurPath ...
    'GLOBCURRENT_geo_MoorSPICE_20120101-20141231_L4_v3.mat'],...
    'u_geo','v_geo')
  lon = double(lon);
  lat = double(lat);
  % Interp satellite geostrophic surface velocities onto mooring location
  % and time
  UV_SC.U = squeeze(interp3(lon,lat,time,permute(u_surface,[2 1 3]),...
    params.lon,params.lat,UV.time))';
  UV_SC.V = squeeze(interp3(lon,lat,time,permute(v_surface,[2 1 3]),...
    params.lon,params.lat,UV.time))';
  UV_SC.U_geo = squeeze(interp3(lon,lat,time,permute(u_geo,[2 1 3]),...
    params.lon,params.lat,UV.time))';
  UV_SC.V_geo = squeeze(interp3(lon,lat,time,permute(v_geo,[2 1 3]),...
    params.lon,params.lat,UV.time))';
  UV_SC.U_ekm = UV_SC.U - UV_SC.U_geo;
  UV_SC.V_ekm = UV_SC.V - UV_SC.V_geo;
  % Finish surface geostrophic velocity structure
  UV_SC.lat = params.lat;
  UV_SC.lon = params.lon;
  UV_SC.time = UV.time;
  % Housekeeping
  clear lat lon time u_surface v_surface
  
  
  %% Save data for surface extrapolation calculation
  fname = [driveName dataPath params.channel '/' params.moor_name ...
    '/UV_GlobCurrents.mat'];
  save(fname,'UV_SC')
  
end