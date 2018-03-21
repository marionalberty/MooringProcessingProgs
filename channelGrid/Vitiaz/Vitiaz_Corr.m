% Vitiaz_Corr.m
% Investigate the correlation bewteen moorings of the sub-inertial velocity
% records. Try to determine if there is a relations that can be applied.

clear all; close all; clc
warning off

set(0,'defaultaxesfontsize',12,'defaultaxeslinewidth',0.7,...
  'defaultlinelinewidth',1,'defaultpatchlinewidth',0.7,...
  'defaultFigureColor','white')


%% Set paths and directories
% Personnal paths:
% driveName='/Users/cyrilgermineaud/Documents/MATLAB/';
% addpath([driveName 'Routines_Cyril'])
driveName ='/Users/marionsofiaalberty/MATLAB/Solomon_Sea/';
dataPath = 'Moorings/Data/Gridded/Vitiaz/';


%% Set parameters
dt = 1;     %[days]
t_filt = 7; %[days]


%% Load mooring data
% Vittiaz East
load([driveName dataPath 'VitiazEast/VitiazEast_linearInterp.mat'])
[UV_E.asv,~] = uvrot(UV.intermediate.V,UV.intermediate.U,-55);
UV_E.time = UV.time;
UV_E.z = UV.intermediate.pres;
params_E = params;
% Vittiaz West
load([driveName dataPath 'VitiazWest/VitiazWest_linearInterp.mat'])
[UV_W.asv,~] = uvrot(UV.intermediate.V,UV.intermediate.U,-55);
UV_W.time = UV.time;
UV_W.z = UV.intermediate.pres;
params_W = params;

% Vittiaz Middle
load([driveName dataPath 'VitiazMiddle/VitiazMiddle_linearInterp.mat'])
UV_M.time = UV.time(1:end-3);
params_M = params;
% Extract Middle velocity at appropriate depths
z_barE = nanmean(UV_E.z);
[~,iz_ME] = min(abs(UV.z-z_barE));
z_barW = nanmean(UV_W.z);
[~,iz_MW] = min(abs(UV.z-z_barW));
[UV_M.asvE,~] = uvrot(UV.V(iz_ME,1:end-3),UV.U(iz_ME,1:end-3),-55);
[UV_M.asvW,~] = uvrot(UV.V(iz_MW,1:end-3),UV.U(iz_MW,1:end-3),-55);
UV_M.zE = UV.z(iz_ME);
UV_M.zW = UV.z(iz_MW);

clear UV params


%% Filter and subsample data to sub-inertial
% Create daily time vector
time = ceil(UV_M.time(1)):floor(UV_M.time(end));

% Filter east and west data
dt_inst = mean(diff(UV_M.time));
w = hamming(round(t_filt/dt_inst));
UV_E.asv = filt_ends(w,UV_E.asv)';
UV_M.asvE = filt_ends(w,UV_M.asvE)';
UV_W.asv = filt_ends(w,UV_W.asv)';
UV_M.asvW = filt_ends(w,UV_M.asvW)';

UV_E.z = filt_ends(w,UV_E.z)';
UV_W.z = filt_ends(w,UV_W.z)';

% Subsample
UV_E.asv = interp1(UV_E.time,UV_E.asv,time);
UV_W.asv = interp1(UV_W.time,UV_W.asv,time);
UV_M.asvE = interp1(UV_M.time,UV_M.asvE,time);
UV_M.asvW = interp1(UV_M.time,UV_M.asvW,time);

UV_E.z = interp1(UV_E.time,UV_E.z,time);
UV_W.z = interp1(UV_W.time,UV_W.z,time);


%% Compare East and Middle

[rW,lagW] = xcorr(UV_W.asv,UV_M.asvW);
[rE,lagE] = xcorr(UV_E.asv,UV_M.asvE);
plot(lagW,rW,'r',lagE,rE,'b')

corrcoef([UV_E.asv' UV_M.asvE'])
corrcoef([UV_W.asv' UV_M.asvW'])

pW = polyfit(UV_M.asvW,UV_W.asv,1);
pE = polyfit(UV_M.asvE,UV_E.asv,1);