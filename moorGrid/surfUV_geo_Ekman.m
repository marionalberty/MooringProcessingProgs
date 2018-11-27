% surfUV_geo_Ekman.m

% Load, interpolate, and save surface geostrophic and Ekman velocities for
% the Solomon Strait moorings.

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
windPath = 'CCMP/MoorSPICE/';
slaPath = 'AVISO/MoorSPICE/';


%% Set parameters
% Coefficients for eddy viscosity empirical formulation determined by 
% Santiago-Mandujano and Firing (1990)
W1 = 1;         %[m/s]
alpha = 8e-5;   %[m^2/s]
b = 2.2;        % exponent
% Density for kinematic wind stress
rho0 = 1019;     %[kg/m^3]


% Choose cases
for i_moor = 1:4
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
    case 4
      % Vitiaz Middle
      params.moor_name = 'VitiazMiddle';
      params.channel = 'Vitiaz';
      load([driveName dataPath params.channel '/' params.moor_name ...
        '/VitiazMiddle.mat'])
      params.lat = params.lat + 0.1;
  end
  
  
  %% Load CCMP winds, interpolate, calc windstress
  load([driveName windPath ...
    'CCMP_MoorSPICE_20120101-20141231_V02.0_L3.0_RSS.mat'])
  lon = double(lon);
  lat = double(lat);
  % Interp wind onto mooring location and time
  uwnd = squeeze(interp3(lon,lat,time,permute(uwnd,[2 1 3]),...
    params.lon,params.lat,UV.time))';
  vwnd = squeeze(interp3(lon,lat,time,permute(vwnd,[2 1 3]),...
    params.lon,params.lat,UV.time))';
  % Calc wind speed
  W = sqrt(uwnd.^2 + vwnd.^2);
  % Calculate eddy viscosity Santiago-Mandujano and Firing (1990)
  A = alpha*(W/W1).^b;
  A(W < 1) = alpha;
  % Calc kinematic wind stress
  [tau_x,tau_y] = ra_windstr(uwnd,vwnd);
  tau_x = tau_x./rho0;
  tau_y = tau_y./rho0;
  % Housekeeping
  clear lon lat time uwnd vwnd
  
  
  %% Calculate the Ekman depth
  % Get Coriolis frequency
  f = sw_f(params.lat);
  % Calc Ekman depth
  d_E = double(sqrt(2*A/abs(f)));     %[m]
  
  
  %% Estimate the Ekman velocity profile for each timestep
  % Make matrix z and tau
  z_E = repmat(UV.z,1,numel(UV.time));
  d_E = repmat(d_E,numel(UV.z),1);
  tau_x = repmat(tau_x,numel(UV.z),1);
  tau_y = repmat(tau_y,numel(UV.z),1);
  % Calculate Ekman spirals
  % Vallis
  UV_E.U = (sqrt(2)./(f*d_E)).*exp(-z_E./d_E).*...
    ((tau_x.*cos((-z_E./d_E)-(pi/4)))-(tau_y.*sin((-z_E./d_E)-(pi/4))));
  UV_E.V = (sqrt(2)./(f*d_E)).*exp(-z_E./d_E).*...
    ((tau_x.*sin((-z_E./d_E)-(pi/4)))+(tau_y.*cos((-z_E./d_E)-(pi/4))));
  % Finish Ekman velocity structure
  UV_E.lat = params.lat;
  UV_E.lon = params.lon;
  UV_E.time = UV.time;
  UV_E.z = UV.z;
  
  
  %% Load AVISO geostrophic velocities, interpolate
  load([driveName slaPath ...
    'SEALEVEL_MoorSPICE_20120101-20141231_L4_v3.mat'],'lat','lon',...
    'time','ugeo','vgeo')
  lon = double(lon);
  lat = double(lat);
  % Interp satellite geostrophic surface velocities onto mooring location
  % and time
  UV_SG.U = squeeze(interp3(lon,lat,time,permute(ugeo,[2 1 3]),...
    params.lon,params.lat,UV.time))';
  UV_SG.V = squeeze(interp3(lon,lat,time,permute(vgeo,[2 1 3]),...
    params.lon,params.lat,UV.time))';
  % Finish surface geostrophic velocity structure
  UV_SG.lat = params.lat;
  UV_SG.lon = params.lon;
  UV_SG.time = UV.time;
  % Housekeeping
  clear lat lon time ugeo vgeo
  
  
  %% Save data for surface extrapolation calculation
  fname = [driveName dataPath params.channel '/' params.moor_name ...
    '/UV_Ekman_SurfGeo.mat'];
  save(fname,'UV_E','UV_SG')
  
end