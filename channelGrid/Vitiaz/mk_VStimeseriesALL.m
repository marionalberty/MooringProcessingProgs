% mk_VStimeseriesALL.m

% Generate filtered, daily timeseries for Vitiaz Strait using the
% different mooring interpolation and extrapolation products.

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
dz = 20;    %[m]
dt = 1;     %[days]
tfilt = 7; %[days]


%% Load bathy and initialize channel grid
% bathy
load([driveName dataPath 'Xsection_Bathy.mat'])

% z-grid
z = [10:dz:max(bathy.z)+dz]';

for i_method = 1:2
  for i_surf = 1:2
    for i_bot = 1:2
      for i_xpass = 1:3
        
        %% Set cases of mooring interpolations to use
        % Mooring interpolation method
        switch i_method
          case 1
            method = 'linearInterp';
          case 2
            method = 'pchipInterp';
        end
        
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
            % Cubic Spline
            xpass = 'pchip';
        end
        
        
        %% Load mooring data
        % Vittiaz East
        load([driveName dataPath 'VitiazEast/VitiazEast_' method '.mat'])
        params.intermediateUV = UV.intermediate;    params_E = params;
        pE = 0.5145;
                
        % Vittiaz Middle
        load([driveName dataPath 'VitiazMiddle/VitiazMiddle_' method...
          '_' bot '_' surf '.mat'])
        UV_M = UV;    SS_M = SS;      TT_M = TT;    params_M = params;
        
        % Vittiaz West
        load([driveName dataPath 'VitiazWest/VitiazWest_' method '.mat'])
        params.intermediateUV = UV.intermediate;    params_W = params;
        pW = 0.6284;
        
        clear UV params SS TT
        
        
        %% Set up grids
        % Create daily time vector
        % Truncate early since up adcp died
        time = ceil(UV_M.time(1)):datenum(2014,2,14);
        
        % velocity grids
        ASV = zeros(numel(z),numel(bathy.dist));
        for i = 1:numel(bathy.dist)
          ASV(z >= bathy.z(i),i) = nan;
        end
        
        % Convert ASV and XSV into 3-D matricies
        ASV = repmat(ASV,1,1,numel(time));
        XSV = ASV;
        
        
        %% Filter, sub-sample, and insert mooring data into channel grid
        % Vitiaz Middle
        % Extract data to be filtered
        UV_filt.U = UV_M.U;       UV_filt.V = UV_M.V;
        % Filter
        UV_filt = moorFilter(UV_filt,UV_M.time,tfilt);
        % Subsample
        U_M = interp1(UV_M.time,UV_filt.U',time)';
        V_M = interp1(UV_M.time,UV_filt.V',time)';
        
        clear UV_filt
        
        
        %% Filer and subsample density data
        % Solomon M1
        % Extract data to be filtered
        ist = 2:9007;
        M.temp = naninterp(TT_M.temp(:,ist)')';
        M.psal = naninterp(SS_M.psal(:,ist)')';
        % Filter
        M_filt = moorFilter(M,TT_M.time(ist),tfilt);
        % Extract full coverage profiles
%         i_s = find(sum(isnan(M_filt.psal),2) == 0);
        i_s = find(sum(isnan(SS_M.psal(:,ist)),2)/...
          numel(SS_M.time(:,ist)) < 0.3);
        % Subsample
        t_M = interp1(TT_M.time(ist),M_filt.temp(i_s,:)',time)';
        s_M = interp1(TT_M.time(ist),M_filt.psal(i_s,:)',time)';
        z_M = TT_M.z(i_s);
        
        clear M* i_s TT_* SS_*
        
        
        %% Calculate density at mooring
        rho_M = sw_dens(s_M,t_M,z_M);
        
        
        %% Get distance between moorings and grid points
        dist_E = nan(size(bathy.z));
        dist_M = nan(size(bathy.z));
        dist_W = nan(size(bathy.z));
        for i = 1:numel(dist_M)
          dist_E(i) = m_lldist([bathy.lon(i) params_E.lon],...
            [bathy.lat(i) params_E.lat]);
          dist_M(i) = m_lldist([bathy.lon(i) params_M.lon],...
            [bathy.lat(i) params_M.lat]);
          dist_W(i) = m_lldist([bathy.lon(i) params_W.lon],...
            [bathy.lat(i) params_W.lat]);
        end
        % Get index of closest grid point to each mooring
        [~,i_dE] = min(dist_E);
        [~,i_dM] = min(dist_M);
        [~,i_dW] = min(dist_W);
        
        % Get angle for ASV & XSV rotation
        phi = azimuth(bathy.lat(1),bathy.lon(1),...
          bathy.lat(end),bathy.lon(end))-90;
        
        %% Put data in grid
        % Middle
        [ASV(1:numel(UV_M.z),i_dM,:),XSV(1:numel(UV_M.z),i_dM,:)] = ...
          uvrot(V_M,U_M,phi);
        % East
        i_zE = find(isnan(ASV(:,i_dE,1)),1,'first')-2;
        ASV(1:i_zE,i_dE,:) = pE*ASV(1:i_zE,i_dM,:);
        % West
        i_zW = find(isnan(ASV(:,i_dW,1)),1,'first')-2;
        ASV(1:i_zW,i_dW,:) = pW*ASV(1:i_zW,i_dM,:);
        
        
        %% Insert t,s,& rho data into full matrix
        RHO = nan(numel(z),numel(bathy.dist),numel(time));
        TEM = nan(numel(z),numel(bathy.dist),numel(time));
        SAL = nan(numel(z),numel(bathy.dist),numel(time));
        
        % Make first column the same as VM
        RHO(z >= z_M(1) & z <= z_M(end),1,:) = rho_M;
        TEM(z >= z_M(1) & z <= z_M(end),1,:) = t_M;
        SAL(z >= z_M(1) & z <= z_M(end),1,:) = s_M;
        
        % Insert VM
        RHO(z >= z_M(1) & z <= z_M(end),i_dM,:) = rho_M;
        TEM(z >= z_M(1) & z <= z_M(end),i_dM,:) = t_M;
        SAL(z >= z_M(1) & z <= z_M(end),i_dM,:) = s_M;
        
        % Make last column the same as VM
        RHO(z >= z_M(1) & z <= z_M(end),end,:) = rho_M;
        TEM(z >= z_M(1) & z <= z_M(end),end,:) = t_M;
        SAL(z >= z_M(1) & z <= z_M(end),end,:) = s_M;
        
        % Interp across channel and add nan sidewalls in
        for k = 1:375
          % Interp across channel
          RHO(:,:,k) = naninterp(RHO(:,:,k)')';
          TEM(:,:,k) = naninterp(TEM(:,:,k)')';
          SAL(:,:,k) = naninterp(SAL(:,:,k)')';
          % Put nans in for bottom boundaries
          for i = 1:numel(bathy.dist)
            RHO(z >= bathy.z(i),i,k) = nan;
            TEM(z >= bathy.z(i),i,k) = nan;
            SAL(z >= bathy.z(i),i,k) = nan;
          end
        end
        
        
        %% Interpolate between moorings and extrapolate to side walls
        [ASV,XSV] = channelInterpVS(ASV,XSV,bathy.dist,[i_dW i_dM i_dE],...
          xpass);
        
        
        %% Plot data
        % save xpassage scheme
        params_M.xpassage_method = ['corrEW_' xpass];
        
        % Collect instrument mean locations
        inst_xz = [bathy.dist(i_dW)*...
          ones(size(nanmean(params_W.intermediateUV.pres,2))) ...
          nanmean(params_W.intermediateUV.pres,2);...
          bathy.dist(i_dM)*...
          ones(size(nanmean(params_M.intermediateUV.pres,2))) ...
          nanmean(params_M.intermediateUV.pres,2);...
          bathy.dist(i_dE)*...
          ones(size(nanmean(params_E.intermediateUV.pres,2))) ...
          nanmean(params_E.intermediateUV.pres,2)];
        
        % Plot
        channelPlot(ASV,XSV,z,time,bathy,params_M,inst_xz,RHO,TEM,SAL)
        close all
        
        
        %% Save data
        f_out = [driveName dataPath params_M.channel '_corrEW_' xpass...
          'Xsection_from_' method '_' bot '_' surf '.mat'];
        save(f_out,'ASV','z','time','bathy','RHO','TEM','SAL','inst_xz')
        
        
      end
    end
  end
end