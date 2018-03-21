% test_STGcompassError.m

clear all
% close all
clc

set(0,'defaultaxesfontsize',12,'defaultaxeslinewidth',0.7,...
  'defaultlinelinewidth',1,'defaultpatchlinewidth',0.7,...
  'defaultFigureColor','white')

% Set case
i_moor = 2;

%% Establish parameters and time
driveName='/Users/marionsofiaalberty/MATLAB/Solomon_Sea/';

params.pathin =[driveName 'Moorings/Data/Raw/'];
params.pathfig = [driveName 'Moorings/Figures/Raw/'];

switch i_moor
  case 1
    % St George's East
    params.moor_name = 'StGeorgesEast';
    params.startgrid_mat = datenum(2012,07,20,08,00,00);
    params.startgrid = datestr(params.startgrid_mat);
    params.stopgrid_mat = datenum(2014,03,08);
    params.stopgrid = datestr(params.stopgrid_mat);
    params.gridDepth = 1430;
  case 2
    % St George's West
    params.moor_name = 'StGeorgesWest';
    params.startgrid_mat = datenum(2012,07,20,01,00,00);
    params.startgrid = datestr(params.startgrid_mat);
    params.stopgrid_mat = datenum(2014,03,07,21,00,00);
    params.stopgrid = datestr(params.stopgrid_mat);
    params.gridDepth = 1240;
end

dt = 1/24;  % [days]
time = params.startgrid_mat:dt:params.stopgrid_mat;


%% Load Up Velocity

% UP
switch i_moor
  case 1
    load([params.pathin params.moor_name ...
      '/ADCP/StGeorgesEast_152m_RDI300kHz_16832_wCE.mat'])
  case 2
    load([params.pathin params.moor_name ...
      '/ADCP/StGeorgesWest_152m_RDI300kHz_16833_wCE.mat'])
end
heading_UP = adcp.heading;
pitch_UP = adcp.pitch;
roll_UP = adcp.roll;
time_UP = adcp.mtime;
nbin_UP = adcp.config.n_cells;
adcp_UP = adcp;
config_UP = config;

% Save copy for later
u_ui = interp1(time_UP,adcp.east_vel',time)';
v_ui = interp1(time_UP,adcp.north_vel',time)';

%% Step through time and rotate all velocity into beam coordinates

% Initialize beam velocity matricies
beam1_vel = nan(size(adcp.east_vel));
beam2_vel = nan(size(adcp.east_vel));
beam3_vel = nan(size(adcp.east_vel));

% loop through time
for i = 1:length(time_UP)
  % Calc inverse of transformation matrix
  % enu = T*beam
  % beam = inv(T)*enu
  T = transformADCPcoord(heading_UP(i),pitch_UP(i),roll_UP(i));
  % loop through each bin at each time step
  for k = 1:nbin_UP
    beam = T\[adcp.east_vel(k,i); adcp.north_vel(k,i); adcp.vert_vel(k,i)];
    beam1_vel(k,i) = beam(1);
    beam2_vel(k,i) = beam(2);
    beam3_vel(k,i) = beam(3);
  end
end


%% Download DOWN data

switch i_moor
  case 1
    % St George's East
    load([params.pathin params.moor_name ...
      '/ADCP/StGeorgesEast_154m_RDI75kHz_16768.mat'])
  case 2
    % St George's West
    load([params.pathin params.moor_name ...
      '/ADCP/StGeorgesWest_154m_RDI75kHz_8866.mat'])
end
heading_DOWN = interp1(adcp.mtime,adcp.heading,time_UP);
pitch_DOWN = interp1(adcp.mtime,adcp.pitch,time_UP);
roll_DOWN = interp1(adcp.mtime,adcp.roll,time_UP);

%% Calculate heading offset from pitch and roll

switch i_moor
  case 1
    % St George's East
    % Use all since no clear issues with pitch and roll
    hoff = fminsearch('checktilt',0,[],[roll_UP; pitch_UP; roll_DOWN; ...
      pitch_DOWN]);
  case 2
    % St George's West
    % Only use a subset since pitch maxes out for UP in the first part of
    % the deployment.
    hoff = fminsearch('checktilt',0,[],[roll_UP(9000:13500); ...
      pitch_UP(9000:13500); roll_DOWN(9000:13500); pitch_DOWN(9000:13500)]);
end

%% Step through time and rotate all velocity back to earth with a compass
%  correction added to the heading

% Initialize beam velocity matricies
east_vel_C = nan(size(beam1_vel));
north_vel_C = nan(size(beam1_vel));
vert_vel_C = nan(size(beam1_vel));

% loop through time
for i = 1:length(time_UP)
  % Calc inverse of transformation matrix
  % enu = T*beam
  % beam = inv(T)*enu
  T = transformADCPcoord(heading_DOWN(i)+hoff,pitch_UP(i),roll_UP(i));
  % loop through each bin at each time step
  for k = 1:nbin_UP
    enu = T * [beam1_vel(k,i); beam2_vel(k,i); beam3_vel(k,i)];
    east_vel_C(k,i) = enu(1);
    north_vel_C(k,i) = enu(2);
    vert_vel_C(k,i) = enu(3);
  end
end

%% Interp new UP data onto regular time grid

dataU.u = interp1(time_UP,east_vel_C',time)';
dataU.v = interp1(time_UP,north_vel_C',time)';

%% Interp DOWN data onto regular time grid

dataD.u = interp1(adcp.mtime,adcp.east_vel',time)';
dataD.v = interp1(adcp.mtime,adcp.north_vel',time)';

%% Concatinate

% Corrected version
U = [flipud(dataU.u); dataD.u];
V = [flipud(dataU.v); dataD.v];
vel = U + 1i*V;
mag = abs(vel);
ang = angle(vel);

% Initial version
U_i = [flipud(u_ui); dataD.u];
V_i = [flipud(v_ui); dataD.v];
vel_i = U_i + 1i*V_i;
mag_i = abs(vel_i);
ang_i = angle(vel_i);


%% Direction Comp plot

figure('position',[0 0 800 400])
ax(1) = subplot(211);
pcolor(ang_i)
shading flat
axis ij
colormap(redblue(25))
colorbar
ylabel('Bin Number')
title([params.moor_name ...
  ' Raw Flow Direction [rad] without compass correction'])
ax(2) = subplot(212);
pcolor(ang)
shading flat
axis ij
colormap(redblue(25))
colorbar
ylabel('Bin Number')
title([params.moor_name ...
  ' Raw Flow Direction [rad] with compass correction'])
xlabel('Time Count')
linkaxes(ax,'xy')
% print figure
figname = [params.pathfig params.moor_name ...
  '/upADCPcompassCorrection_direction.png'];
print(figname,'-dpng')


%% Magnitude Comp plot

figure('position',[0 0 800 400])
ax(1) = subplot(211);
pcolor(mag_i)
shading flat
axis ij
colormap(parula(25))
colorbar
ylabel('Bin Number')
title([params.moor_name ...
  ' Raw Flow Speed [m/s] without compass correction'])
ax(2) = subplot(212);
pcolor(mag)
shading flat
axis ij
colormap(parula(25))
colorbar
ylabel('Bin Number')
title([params.moor_name ...
  ' Raw Flow Speed [m/s] with compass correction'])
xlabel('Time Count')
linkaxes(ax,'xy')
% print figure
figname = [params.pathfig params.moor_name ...
  '/upADCPcompassCorrection_magnitude.png'];
print(figname,'-dpng')