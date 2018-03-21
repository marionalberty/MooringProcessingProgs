% mk_SStimeseries.m

% Generate filtered, daily timeseries for Solomon Strait using the
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
dataPath = 'Moorings/Data/Gridded/SolomonStrait/';


%% Set parameters
dz = 20;    %[m]
dt = 1;     %[days]
tfilt = 7; %[days]


%% Load bathy and initialize channel grid
% bathy
load([driveName dataPath 'Xsection_Bathy.mat'])

% z-grid
z = [10:dz:max(bathy.z)]';


%% Do interpolation with velocity
for i_method = 1:2
  for i_surf = 1:2
    for i_bot = 1:2
      for i_xpass = 1:2
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
        end
        
        
        %% Load mooring data
        % Solomon M1
        load([driveName dataPath 'Solomon_M1/Solomon_M1_' method ...
          '_' bot '_' surf '.mat'])
        UV_1 = UV;    SS_1 = SS;      TT_1 = TT;    params_1 = params;
        
        % Solomon M2b
        load([driveName dataPath 'Solomon_M2b/Solomon_M2b_' method ...
          '_' bot '_' surf '.mat'])
        UV_2 = UV;    SS_2 = SS;      TT_2 = TT;    params_2 = params;
        
        % Solomon M3
        load([driveName dataPath 'Solomon_M3/Solomon_M3_' method ...
          '_' bot '_' surf '.mat'])
        UV_3 = UV;    SS_3 = SS;      TT_3 = TT;    params_3 = params;
        
        clear UV params SS TT
        
        
        %% Set up grids
        % Create daily time vector
        time = ceil(max([UV_1.time(1) UV_2.time(1) UV_3.time(1)])):...
          floor(min([UV_1.time(end) UV_2.time(end) UV_3.time(end)]));
        
        % velocity grids
        ASV = zeros(numel(z),numel(bathy.dist));
        for i = 1:numel(bathy.dist)
          ASV(z >= bathy.z(i),i) = nan;
        end
        
        % Convert ASV and XSV into 3-D matricies
        ASV = repmat(ASV,1,1,numel(time));
        XSV = ASV;
        
        
        %% Filter, sub-sample, and insert velocity data into channel grid
        % Solomon M1
        % Extract data to be filtered
        UV_filt.U = UV_1.U;       UV_filt.V = UV_1.V;
        % Filter
        UV_filt = moorFilter(UV_filt,UV_1.time,tfilt);
        % Subsample
        U_1 = interp1(UV_1.time,UV_filt.U',time)';
        V_1 = interp1(UV_1.time,UV_filt.V',time)';
        
        % Solomon M2
        % Extract data to be filtered
        UV_filt.U = UV_2.U;       UV_filt.V = UV_2.V;
        % Filter
        UV_filt = moorFilter(UV_filt,UV_2.time,tfilt);
        % Subsample
        U_2 = interp1(UV_2.time,UV_filt.U',time)';
        V_2 = interp1(UV_2.time,UV_filt.V',time)';
        
        % Solomon M3
        % Extract data to be filtered
        UV_filt.U = UV_3.U;       UV_filt.V = UV_3.V;
        % Filter
        UV_filt = moorFilter(UV_filt,UV_3.time,tfilt);
        % Subsample
        U_3 = interp1(UV_3.time,UV_filt.U',time)';
        V_3 = interp1(UV_3.time,UV_filt.V',time)';
        
        clear UV_filt
        
        
        %% Filer and subsample density data
        % Solomon M1
        % Extract data to be filtered
        M1.temp = naninterp(TT_1.temp')';   M1.psal = naninterp(SS_1.psal')';
        % Filter
        M_filt = moorFilter(M1,TT_1.time,tfilt);
        % Extract full coverage profiles
        i_s = find(sum(isnan(M_filt.psal),2) == 0);
        % Subsample
        t_1 = interp1(TT_1.time,M_filt.temp(i_s,:)',time)';
        s_1 = interp1(TT_1.time,M_filt.psal(i_s,:)',time)';
        z_1 = TT_1.z(i_s);
        
        % Solomon M2
        % Extract data to be filtered
        M2.temp = naninterp(TT_2.temp')';   M2.psal = naninterp(SS_2.psal')';
        % Filter
        M_filt = moorFilter(M2,TT_2.time,tfilt);
        % Extract full coverage profiles
        i_s = find(sum(isnan(M_filt.psal),2) == 0);
        % Subsample
        t_2 = interp1(TT_2.time,M_filt.temp(i_s,:)',time)';
        s_2 = interp1(TT_2.time,M_filt.psal(i_s,:)',time)';
        z_2 = TT_2.z(i_s);
        
        % Solomon M3
        % Extract data to be filtered
        M3.temp = naninterp(TT_3.temp')';   M3.psal = naninterp(SS_3.psal')';
        % Filter
        M_filt = moorFilter(M3,TT_3.time,tfilt);
        % Extract full coverage profiles
        i_s = find(sum(isnan(M_filt.psal),2) == 0);
        % Subsample
        t_3 = interp1(TT_3.time,M_filt.temp(i_s,:)',time)';
        s_3 = interp1(TT_3.time,M_filt.psal(i_s,:)',time)';
        z_3 = TT_3.z(i_s);
        
        clear M* i_s TT_* SS_*
        
        
        %% Calculate density at each mooring
        rho_1 = sw_dens(s_1,t_1,0);
        rho_2 = sw_dens(s_2,t_2,0);
        rho_3 = sw_dens(s_3,t_3,0);
        
        
        %% Get the depth of the 27 & 27.5 sig theta isopycnal
        z1_sgth267 = nan(size(time));
        z2_sgth267 = nan(size(time));
        z3_sgth267 = nan(size(time));
        z1_sgth275 = nan(size(time));
        z2_sgth275 = nan(size(time));
        z3_sgth275 = nan(size(time));
        % Cycle through time
        for i = 1:numel(time)
          z1_sgth267(i) = interp1(rho_1(:,i),z_1,1026.7);
          z2_sgth267(i) = interp1(rho_2(:,i),z_2,1026.7);
          z3_sgth267(i) = interp1(rho_3(:,i),z_3,1026.7);
          z1_sgth275(i) = interp1(rho_1(:,i),z_1,1027.5);
          z2_sgth275(i) = interp1(rho_2(:,i),z_2,1027.5);
          z3_sgth275(i) = interp1(rho_3(:,i),z_3,1027.5);
        end
        
        
        %% Get mooring locations and segment angles
        % Get distance between moorings and grid points
        dist_1 = nan(size(bathy.z));
        dist_2 = nan(size(bathy.z));
        dist_3 = nan(size(bathy.z));
        for i = 1:numel(dist_1)
          dist_1(i) = m_lldist([bathy.lon(i) params_1.lon],...
            [bathy.lat(i) params_1.lat]);
          dist_2(i) = m_lldist([bathy.lon(i) params_2.lon],...
            [bathy.lat(i) params_2.lat]);
          dist_3(i) = m_lldist([bathy.lon(i) params_3.lon],...
            [bathy.lat(i) params_3.lat]);
        end
        % Get index of closest grid point to each mooring
        [~,i_d1] = min(dist_1);
        [~,i_d2] = min(dist_2);
        [~,i_d3] = min(dist_3);
        
        % Get angles for ASV & XSV rotation
        phi_12 = azimuth(bathy.lat(1),bathy.lon(1),...
          bathy.lat(i_d2),bathy.lon(i_d2))-90;
        phi_23 = azimuth(bathy.lat(i_d2),bathy.lon(i_d2),...
          bathy.lat(end),bathy.lon(end))-90;
        
        clear dist_*
        
        
        %% Put velocity data in grid
        % M1
        [ASV(1:numel(UV_1.z),i_d1,:),XSV(1:numel(UV_1.z),i_d1,:)] = ...
          uvrot(V_1,U_1,phi_12);
        
        % M2
        [ASV(1:numel(UV_2.z),i_d2-1,:),XSV(1:numel(UV_2.z),i_d2-1,:)] = ...
          uvrot(V_2,U_2,phi_12);
        [ASV(1:numel(UV_2.z),i_d2,:),XSV(1:numel(UV_2.z),i_d2,:)] = ...
          uvrot(V_2,U_2,phi_23);
        
        % M3
        [ASV(1:numel(UV_3.z),i_d3,:),XSV(1:numel(UV_3.z),i_d3,:)] = ...
          uvrot(V_3,U_3,phi_23);
        
        
        %% Insert t,s,& rho data into full matrix
        RHO = nan(numel(z),numel(bathy.dist),numel(time));
        TEM = nan(numel(z),numel(bathy.dist),numel(time));
        SAL = nan(numel(z),numel(bathy.dist),numel(time));
        
        % Make first column the same as M1
        RHO(z >= z_1(1) & z <= z_1(end),1,:) = rho_1;
        TEM(z >= z_1(1) & z <= z_1(end),1,:) = t_1;
        SAL(z >= z_1(1) & z <= z_1(end),1,:) = s_1;
        
        % Insert M1
        RHO(z >= z_1(1) & z <= z_1(end),i_d1,:) = rho_1;
        TEM(z >= z_1(1) & z <= z_1(end),i_d1,:) = t_1;
        SAL(z >= z_1(1) & z <= z_1(end),i_d1,:) = s_1;
        
        % Insert M2
        RHO(z >= z_2(1) & z <= z_2(end),i_d2,:) = rho_2;
        TEM(z >= z_2(1) & z <= z_2(end),i_d2,:) = t_2;
        SAL(z >= z_2(1) & z <= z_2(end),i_d2,:) = s_2;
        
        % Insert M3
        RHO(z >= z_3(1) & z <= z_3(end),i_d3,:) = rho_3;
        TEM(z >= z_3(1) & z <= z_3(end),i_d3,:) = t_3;
        SAL(z >= z_3(1) & z <= z_3(end),i_d3,:) = s_3;
        
        % Make last column the same as M3
        RHO(z >= z_3(1) & z <= z_3(end),end,:) = rho_3;
        TEM(z >= z_3(1) & z <= z_3(end),end,:) = t_3;
        SAL(z >= z_3(1) & z <= z_3(end),end,:) = s_3;
        
        % Interp across channel and add nan sidewalls in
        for k = 1:numel(time)
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
        [ASV,XSV] = channelInterpSS(ASV,XSV,bathy.dist,[i_d1 i_d2-1 ...
          i_d2 i_d3],xpass);
        
        
        %% Interpolate z_sgth27 across channel
        z_sgth267 = [z1_sgth267' z1_sgth267' z2_sgth267' z3_sgth267' z3_sgth267'];
        z_sgth267 = interp1(bathy.dist([1 i_d1 i_d2 i_d3 end]),z_sgth267',...
          bathy.dist,xpass)';
        z_sgth275 = [z1_sgth275' z1_sgth275' z2_sgth275' z3_sgth275' ...
          z3_sgth275'];
        z_sgth275 = interp1(bathy.dist([1 i_d1 i_d2 i_d3 end]),z_sgth275',...
          bathy.dist,xpass)';
        
        
        %% Plot data
        
        % save xpassage scheme
        params_1.xpassage_method = xpass;
        
        % Collect instrument mean locations
        inst_xz = [bathy.dist(i_d1)*...
          ones(size(nanmean(params_1.intermediateUV.pres,2))) ...
          nanmean(params_1.intermediateUV.pres,2);...
          bathy.dist(i_d2)*...
          ones(size(nanmean(params_2.intermediateUV.pres,2))) ...
          nanmean(params_2.intermediateUV.pres,2);...
          bathy.dist(i_d3)*...
          ones(size(nanmean(params_3.intermediateUV.pres,2))) ...
          nanmean(params_3.intermediateUV.pres,2)];
        
        % Plot
        channelPlot(ASV,XSV,z,time,bathy,params_1,inst_xz,RHO,TEM,SAL)
        close all
        
        
        %% Save data
        f_out = [driveName dataPath params_1.channel '_' xpass 'Xsection_'...
          'from_' method '_' bot '_' surf '.mat'];
        save(f_out,'ASV','XSV','z','time','bathy','z_sgth267','z_sgth275',...
          'RHO','TEM','SAL','inst_xz')
        
        
      end
    end
  end
end