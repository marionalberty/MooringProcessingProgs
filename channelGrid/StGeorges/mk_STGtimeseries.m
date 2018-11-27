% mk_STGtimeseries.m

% Generate filtered, daily timeseries for St. George's Channel using the
% different mooring interpolation and extrapolation products.

clear all; close all; clc
warning off

set(0,'defaultaxesfontsize',16,'defaultaxeslinewidth',0.7,...
  'defaultlinelinewidth',1,'defaultpatchlinewidth',0.7,...
  'defaultFigureColor','white')


%% Set paths and directories
% Personnal paths:
% driveName='/Users/cyrilgermineaud/Documents/MATLAB/';
% addpath([driveName 'Routines_Cyril'])
driveName ='/Users/marionsofiaalberty/MATLAB/Solomon_Sea/';
dataPath = 'Moorings/Data/Gridded/StGeorges/';


%% Set parameters
dz = 20;    %[m]
dt = 1;     %[days]
tfilt = 7; %[days]


%% Load bathy and initialize channel grid
% bathy
load([driveName dataPath 'Xsection_Bathy.mat'])

% z-grid
z = transpose(10:dz:max(bathy.z));

for i_surf = 1:2
  for i_bot = 1:2
    for i_xpass = 1:3
      %% Set cases of mooring interpolations to use
      % Mooring surface extrapolation method
      switch i_surf
        case 1
          surf = 'linear2surf';
        case 2
          surf = 'slab2surf';
      end
      
      % Mooring bottom extrapolation method
      switch i_bot
        case 1
          bot = 'linear2bot';
        case 2
          bot = 'slab2bot';
      end
      
      % Cross-passage interpolation scheme
      switch i_xpass
        case 1
          % Linear across passage
          xpass = 'linear';
        case 2
          % Nearest
          xpass = 'nearest';
        case 3
          % Western
          xpass = 'western';
      end
      
      
      %% Load mooring data
      % St George's East
      load([driveName dataPath 'StGeorgesEast/StGeorgesEast_' bot '_' ...
        surf '.mat'])
      UV_E = UV;    SS_E = SS;      TT_E = T4S;      RHO_E = RHO;
      params_E = params;
      
      % St George's West
      load([driveName dataPath 'StGeorgesWest/StGeorgesWest_' bot '_' ...
        surf '.mat'])
      UV_W = UV;    SS_W = SS;      TT_W = T4S;      RHO_W = RHO;
      params_W = params;
      
      % Housekeeping
      clear UV params SS TT RHO T4S
      
      
      %% Set up velocity grids
      % Create daily time vector
      time = ceil(max([UV_E.time(1) UV_W.time(1)])):...
        floor(min([UV_E.time(end) UV_W.time(end)]));
      
      % velocity grids
      ASV = zeros(numel(z),numel(bathy.dist));
      for i = 1:numel(bathy.dist)
        ASV(z >= bathy.z(i),i) = nan;
      end
      
      % Convert ASV and XSV into 3-D matricies
      ASV = repmat(ASV,1,1,numel(time));
      XSV = ASV;
      
      
      %% Filter, sub-sample, and insert velocity data into channel grid
      % St George's East
      % Extract data to be filtered
      UV_filt.U = UV_E.U;         UV_filt.V = UV_E.V;
      % Filter
      UV_filt = moorFilter(UV_filt,UV_E.time,tfilt);
      % Subsample
      U_E = interp1(UV_E.time,UV_filt.U',time)';
      V_E = interp1(UV_E.time,UV_filt.V',time)';
      
      % St George's West
      % Extract data to be filtered
      UV_filt.U = UV_W.U;         UV_filt.V = UV_W.V;
      % Filter
      UV_filt = moorFilter(UV_filt,UV_W.time,tfilt);
      % Subsample
      U_W = interp1(UV_W.time,UV_filt.U',time)';
      V_W = interp1(UV_W.time,UV_filt.V',time)';
      
      % Housekeeping
      clear UV_filt
      
      
      %% Filer and subsample density data
      % St George's East
      ist = 1:9216;
      % Extract data to be filtered
      ME.temp = naninterp(TT_E.temp(:,ist)')';
      ME.psal = naninterp(SS_E.psal(:,ist)')';
      ME.sgth = naninterp(RHO_E.sgth(:,ist)')';
      % Filter
      M_filt = moorFilter(ME,TT_E.time(ist),tfilt);
      % Extract full coverage profiles
      %         i_s = find(sum(isnan(M_filt.psal),2) == 0);
      i_s = find(sum(isnan(SS_E.psal(:,ist)),2)/...
        numel(SS_E.time(:,ist)) < 0.3);
      % Subsample
      t_E = interp1(TT_E.time(ist),M_filt.temp(i_s,:)',time)';
      s_E = interp1(TT_E.time(ist),M_filt.psal(i_s,:)',time)';
      rho_E = interp1(TT_E.time(ist),M_filt.sgth(i_s,:)',time)';
      z_E = TT_E.z(i_s);
      
      % St George's West
      ist = 1:5827;
      % Extract data to be filtered
      MW.temp = naninterp(TT_W.temp(:,ist)')';
      MW.psal = naninterp(SS_W.psal(:,ist)')';
      MW.sgth = naninterp(RHO_W.sgth(:,ist)')';
      % Filter
      M_filt = moorFilter(MW,TT_W.time(ist),tfilt);
      % Extract full coverage profiles
      %         i_s = find(sum(isnan(M_filt.psal),2) == 0);
      i_s = find(sum(isnan(SS_W.psal(:,ist)),2)/...
        numel(SS_W.time(:,ist)) < 0.3);
      % Subsample
      t_W = interp1(TT_W.time(ist),M_filt.temp(i_s,:)',time)';
      s_W = interp1(TT_W.time(ist),M_filt.psal(i_s,:)',time)';
      rho_W = interp1(TT_W.time(ist),M_filt.sgth(i_s,:)',time)';
      z_W = TT_W.z(i_s);
      
      % Housekeeping
      clear M* i_s TT_* SS_* ist
      
      
      %% Get distance between moorings and grid points
      dist_E = nan(size(bathy.z));
      dist_W = nan(size(bathy.z));
      for i = 1:numel(dist_W)
        dist_E(i) = m_lldist([bathy.lon(i) params_E.lon],...
          [bathy.lat(i) params_E.lat]);
        dist_W(i) = m_lldist([bathy.lon(i) params_W.lon],...
          [bathy.lat(i) params_W.lat]);
      end
      % Get index of closest grid point to each mooring
      [~,i_dE] = min(dist_E);
      [~,i_dW] = min(dist_W);
      
      % Get angle for ASV & XSV rotation
      phi = 90 - azimuth(bathy.lat(1),bathy.lon(1),...
        bathy.lat(end),bathy.lon(end));
      
      
      %% Rotate and put data in grid
      % West
      [XSV(1:numel(UV_W.z),i_dW,:),ASV(1:numel(UV_W.z),i_dW,:)] = ...
        uvrot(U_W,V_W,phi);
      
      % East
      [XSV(1:numel(UV_E.z),i_dE,:),ASV(1:numel(UV_E.z),i_dE,:)] = ...
        uvrot(U_E,V_E,phi);
      
      
      %% Insert t,s,& rho data into full matrix
      RHO = nan(numel(z),numel(bathy.dist),numel(time));
      TEM = nan(numel(z),numel(bathy.dist),numel(time));
      SAL = nan(numel(z),numel(bathy.dist),numel(time));
      
      % Make first column the same as SGW
      RHO(z >= z_W(1) & z <= z_W(end),1,:) = rho_W;
      TEM(z >= z_W(1) & z <= z_W(end),1,:) = t_W;
      SAL(z >= z_W(1) & z <= z_W(end),1,:) = s_W;
      
      % Insert SGW
      RHO(z >= z_W(1) & z <= z_W(end),i_dW,:) = rho_W;
      TEM(z >= z_W(1) & z <= z_W(end),i_dW,:) = t_W;
      SAL(z >= z_W(1) & z <= z_W(end),i_dW,:) = s_W;
      
      % Insert SGE
      RHO(z >= z_E(1) & z <= z_E(end),i_dE,:) = rho_E;
      TEM(z >= z_E(1) & z <= z_E(end),i_dE,:) = t_E;
      SAL(z >= z_E(1) & z <= z_E(end),i_dE,:) = s_E;
      
      % Make last column the same as SGE
      RHO(z >= z_E(1) & z <= z_E(end),end,:) = rho_E;
      TEM(z >= z_E(1) & z <= z_E(end),end,:) = t_E;
      SAL(z >= z_E(1) & z <= z_E(end),end,:) = s_E;
      
      
      %% Interpolate between moorings and extrapolate to side walls
      [ASV,XSV] = channelInterpSTG(ASV,XSV,bathy.dist,[i_dW i_dE],...
        xpass);
      
      
      %% Plot data
      
      % save xpassage scheme
      params_E.xpassage_method = xpass;
      
      % Collect instrument mean locations
      inst_xz = [bathy.dist(i_dE)*...
        ones(size(nanmean(params_E.intermediateUV.pres,2))) ...
        nanmean(params_E.intermediateUV.pres,2);...
        bathy.dist(i_dW)*...
        ones(size(nanmean(params_W.intermediateUV.pres,2))) ...
        nanmean(params_W.intermediateUV.pres,2)];
      
      % Plot
      channelPlot(ASV,XSV,z,time,bathy,params_E,inst_xz,RHO,TEM,SAL)
      close all
      
      
      %% Save data
      i_moor = [i_dW i_dE];
      f_out = [driveName dataPath params_E.channel '_' xpass ...
        'Xsection_' bot '_' surf '.mat'];
      save(f_out,'ASV','XSV','z','time','bathy','RHO','TEM','SAL',...
        'inst_xz','i_moor')
      
      
    end
  end
end