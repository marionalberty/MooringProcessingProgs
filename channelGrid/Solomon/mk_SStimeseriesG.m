% mk_SStimeseriesG.m

% Generate filtered, daily timeseries for Solomon Strait using geostrophy
% to estimate velocities between the moorings reference from the bottom.
% See archived version for top referenced option as well which doesn't work
% as well because of missing shallow density observations.

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
g = -9.8;   %[m/s^2]

% Set case for reference velocity
Gref = 'Bot';


%% Load bathy and initialize channel grid
% bathy
load([driveName dataPath 'Xsection_Bathy.mat'])

% z-grid
z = transpose(10:dz:max(bathy.z));


for i_surf = 1:4
  for i_bot = 1:2
    for i_xpass = 1:3
      %% Set cases of mooring interpolations to use
      % Mooring surface extrapolation method
      switch i_surf
        case 1
          surf = 'linear2surf';
        case 2
          surf = 'slab2surf';
        case 3
          surf = 'sat4surf';
        case 4
          surf = 'gcur4surf';
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
      % Solomon M1
      load([driveName dataPath 'Solomon_M1/Solomon_M1_' bot '_' surf ...
        '.mat'])
      UV_1 = UV;    SS_1 = SS;      TT_1 = T4S;    RHO_1 = RHO;
      params_1 = params;
      
      % Solomon M2a and Solomon M2b
      load([driveName dataPath 'Solomon_M2a/Solomon_M2a.mat'])
      SS_2 = SS;    TT_2 = T4S;     RHO_21 = RHO;
      load([driveName dataPath 'Solomon_M2b/Solomon_M2b_' bot '_' surf ...
        '.mat'])
      UV_2 = UV;    params_2 = params;
      
      % Solomon M3
      load([driveName dataPath 'Solomon_M3/Solomon_M3_' bot '_' surf ...
        '.mat'])
      UV_3 = UV;    SS_3 = SS;      TT_3 = T4S;    RHO_3 = RHO;
      params_3 = params;
      
      % Housekeeping
      clear UV TT SS T4S RHO params
      
      
%       %% Adjust RHO grids to only use sensor pairs
%       [RHO_1,RHO_21] = gridRho4Geo(RHO_1,RHO_2);
%       [RHO_22,RHO_3] = gridRho4Geo(RHO_2,RHO_3);
%       % House keeping
%       clear RHO_2
      
      
      %% Set up grid
      % Create daily time vector
      time = ceil(max([TT_1.time(1) TT_2.time(1) TT_3.time(1)])):...
        floor(min([TT_1.time(end) TT_2.time(end) TT_3.time(end)]));
      
      % velocity grid
      ASV_G = zeros(numel(z),numel(bathy.dist));
      for i = 1:numel(bathy.dist)
        ASV_G(z >= bathy.z(i),i) = nan;
      end
      % Convert ASV into 3-D matrix
      ASV_G = repmat(ASV_G,1,1,numel(time));
      
      
      %% Filer and subsample data
      % Solomon M1
      % Extract data to be filtered
      M1.temp = naninterp(TT_1.temp')';
      M1.psal = naninterp(SS_1.psal')';
      M1.sgth = naninterp(RHO_1.sgth')';
      % Filter
      M_filt = moorFilter(M1,TT_1.time,tfilt);
      U_filt = moorFilter(UV_1,UV_1.time,tfilt);
      % Extract full coverage profiles
      i_s = find(sum(isnan(M_filt.sgth),2) == 0);
      % Subsample
      t_1 = interp1(TT_1.time,M_filt.temp(i_s,:)',time)';
      s_1 = interp1(TT_1.time,M_filt.psal(i_s,:)',time)';
      rho_1 = interp1(TT_1.time,M_filt.sgth(i_s,:)',time)';
      z_1 = TT_1.z(i_s);
      u_1 = interp1(UV_1.time,U_filt.U',time)';
      v_1 = interp1(UV_1.time,U_filt.V',time)';
      switch i_surf
        case 3
          ug_1 = interp1(UV_1.time,U_filt.U_geo',time)';
          vg_1 = interp1(UV_1.time,U_filt.V_geo',time)';
          ue_1 = interp1(UV_1.time,(U_filt.U-U_filt.U_geo)',time)';
          ve_1 = interp1(UV_1.time,(U_filt.V-U_filt.V_geo)',time)';
        case 4
          us_1 = interp1(UV_1.time,U_filt.U(1,:)',time)';
          vs_1 = interp1(UV_1.time,U_filt.V(1,:)',time)';
      end
      zu_1 = UV_1.z;
      
      % Solomon M2
      % Extract data to be filtered
      M2.temp = naninterp(TT_2.temp')';
      M2.psal = naninterp(SS_2.psal')';
      M2.sgth = naninterp(RHO_21.sgth')';
%       M2.sgth2 = naninterp(RHO_22.sgth')';
      % Filter
      M_filt = moorFilter(M2,TT_2.time,tfilt);
      U_filt = moorFilter(UV_2,UV_2.time,tfilt);
      % Extract full coverage profiles
      i_s = find(sum(isnan(M_filt.sgth),2) == 0);
%       i_s2 = find(sum(isnan(M_filt.sgth2),2) == 0);
      % Subsample
      t_2 = interp1(TT_2.time,M_filt.temp(i_s,:)',time)';
      s_2 = interp1(TT_2.time,M_filt.psal(i_s,:)',time)';
      rho_21 = interp1(TT_2.time,M_filt.sgth(i_s,:)',time)';
%       rho_22 = interp1(TT_2.time,M_filt.sgth2(i_s2,:)',time)';
      z_21 = TT_2.z(i_s);
%       z_22 = TT_2.z(i_s2);
      u_2 = interp1(UV_2.time,U_filt.U',time)';
      v_2 = interp1(UV_2.time,U_filt.V',time)';
      switch i_surf
        case 3
          ug_2 = interp1(UV_2.time,U_filt.U_geo',time)';
          vg_2 = interp1(UV_2.time,U_filt.V_geo',time)';
          ue_2 = interp1(UV_2.time,(U_filt.U-U_filt.U_geo)',time)';
          ve_2 = interp1(UV_2.time,(U_filt.V-U_filt.V_geo)',time)';
        case 4
          us_2 = interp1(UV_2.time,U_filt.U(1,:)',time)';
          vs_2 = interp1(UV_2.time,U_filt.V(1,:)',time)';
      end
      zu_2 = UV_2.z;
      
      % Solomon M3
      % Extract data to be filtered
      M3.temp = naninterp(TT_3.temp')';
      M3.psal = naninterp(SS_3.psal')';
      M3.sgth = naninterp(RHO_3.sgth')';
      % Filter
      M_filt = moorFilter(M3,TT_3.time,tfilt);
      U_filt = moorFilter(UV_3,UV_3.time,tfilt);
      % Extract full coverage profiles
      i_s = find(sum(isnan(M_filt.sgth),2) == 0);
      % Subsample
      t_3 = interp1(TT_3.time,M_filt.temp(i_s,:)',time)';
      s_3 = interp1(TT_3.time,M_filt.psal(i_s,:)',time)';
      rho_3 = interp1(TT_3.time,M_filt.sgth(i_s,:)',time)';
      z_3 = TT_3.z(i_s);
      u_3 = interp1(UV_3.time,U_filt.U',time)';
      v_3 = interp1(UV_3.time,U_filt.V',time)';
      switch i_surf
        case 3
          ug_3 = interp1(UV_3.time,U_filt.U_geo',time)';
          vg_3 = interp1(UV_3.time,U_filt.V_geo',time)';
          ue_3 = interp1(UV_3.time,(U_filt.U-U_filt.U_geo)',time)';
          ve_3 = interp1(UV_3.time,(U_filt.V-U_filt.V_geo)',time)';
        case 4
          us_3 = interp1(UV_3.time,U_filt.U(1,:)',time)';
          vs_3 = interp1(UV_3.time,U_filt.V(1,:)',time)';
      end
      zu_3 = UV_3.z;
      
      % Housekeeping
      clear M* i_u i_s i_s2 TT_* UV_*
      
      
      %% Get the depth of the 26.7 & 27.5 sig theta isopycnal
      z1_sgth267 = nan(size(time));
      z2_sgth267 = nan(size(time));
      z3_sgth267 = nan(size(time));
      z1_sgth275 = nan(size(time));
      z2_sgth275 = nan(size(time));
      z3_sgth275 = nan(size(time));
      % Cycle through time
      for i = 1:numel(time)
        z1_sgth267(i) = interp1(rho_1(:,i),z_1,1026.7);
        z2_sgth267(i) = interp1(rho_21(:,i),z_21,1026.7);
        z3_sgth267(i) = interp1(rho_3(:,i),z_3,1026.7);
        z1_sgth275(i) = interp1(rho_1(:,i),z_1,1027.5);
        z2_sgth275(i) = interp1(rho_21(:,i),z_21,1027.5);
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
      
      % Get angles for ASV rotation
      phi_12 = 90 - azimuth(bathy.lat(1),bathy.lon(1),...
        bathy.lat(i_d2),bathy.lon(i_d2));
      phi_23 = 90 - azimuth(bathy.lat(i_d2),bathy.lon(i_d2),...
        bathy.lat(end),bathy.lon(end));
      
      % Housekeeping
      clear dist_*
      
      
      %% Insert t,s,& rho data into full matrix
      RHO = nan(numel(z),numel(bathy.dist),numel(time));
      TEM = nan(numel(z),numel(bathy.dist),numel(time));
      SAL = nan(numel(z),numel(bathy.dist),numel(time));
      
      % Insert M1
      RHO(z >= z_1(1) & z <= z_1(end),i_d1,:) = rho_1;
      TEM(z >= z_1(1) & z <= z_1(end),i_d1,:) = t_1;
      SAL(z >= z_1(1) & z <= z_1(end),i_d1,:) = s_1;
      
      % Insert M2
      RHO(z >= z_21(1) & z <= z_21(end),i_d2-1,:) = rho_21;
      TEM(z >= z_21(1) & z <= z_21(end),i_d2-1,:) = t_2;
      SAL(z >= z_21(1) & z <= z_21(end),i_d2-1,:) = s_2;
      
      % Insert M2
      RHO(z >= z_21(1) & z <= z_21(end),i_d2,:) = rho_21;
      TEM(z >= z_21(1) & z <= z_21(end),i_d2,:) = t_2;
      SAL(z >= z_21(1) & z <= z_21(end),i_d2,:) = s_2;
      
      % Insert M3
      RHO(z >= z_3(1) & z <= z_3(end),i_d3,:) = rho_3;
      TEM(z >= z_3(1) & z <= z_3(end),i_d3,:) = t_3;
      SAL(z >= z_3(1) & z <= z_3(end),i_d3,:) = s_3;
      
      
      %% Calculate vertical profile of shear for each line
      % West Segment (M1 & M2, line 1)
      % drho/dx
      z_top1 = max([z_1(1) z_21(1)]);
      z_bot1 = min([z_1(end) z_21(end) bathy.z(i_d1:i_d2)]);
      i_zt = find(z == z_top1);
      i_zb = find(z <= z_bot1,1,'last');
      z_bot1 = z(i_zb);
      drho1 = rho_21(z_21 >= z_top1 & z_21 <= z_bot1,:) - ...
        rho_1(z_1 >= z_top1 & z_1 <= z_bot1,:);
      dx1 = (bathy.dist(i_d2) - bathy.dist(i_d1))*1000;
      % grf = g/rho_0 f
      grf1 = repmat(g/sw_f(mean(bathy.lat(i_d1:i_d2)))./mean(...
        [rho_21(z_21 >= z_top1 & z_21 <= z_bot1,:); ...
        rho_1(z_1 >= z_top1 & z_1 <= z_bot1,:)]),i_zb-i_zt+1,1);
      % dv
      dv1 = dz * grf1 .* drho1 ./ dx1;
      
      % East Segment (M2 & M3, line 2)
      % drho/dx
      z_top2 = max([z_21(1) z_3(1)]);
      z_bot2 = min([z_21(end) z_3(end)]);
      i_zt = find(z == z_top2);
      i_zb = find(z == z_bot2);
      drho2 = rho_3(z_3 >= z_top2 & z_3 <= z_bot2,:) - ...
        rho_21(z_21 >= z_top2 & z_21 <= z_bot2,:);
      dx2 = (bathy.dist(i_d3) - bathy.dist(i_d2))*1000;
      % grf = g/rho_0 f
      grf2 = repmat(g/sw_f(mean(bathy.lat(i_d2:i_d3)))./mean(...
        [rho_3(z_3 >= z_top2 & z_3 <= z_bot2,:); ...
        rho_21(z_21 >= z_top2 & z_21 <= z_bot2,:)]),i_zb-i_zt+1,1);
      % dv
      dv2 = dz * grf2 .* drho2 ./ dx2;
      
      
      %% Use near bottom velocities for reference
      % M1
      [~,asv] = uvrot(u_1,v_1,phi_12);
      asv11 = [cumsum([dv1; asv(zu_1 == z_bot1+dz,:)],'reverse'); ...
        asv(zu_1 > z_bot1+dz,:)];
      z11 = transpose(z_top1:dz:zu_1(end));
      % M2 line 1
      [~,asv] = uvrot(u_2,v_2,phi_12);
      asv12 = [cumsum([dv1; asv(zu_2 == z_bot1+dz,:)],'reverse'); ...
        asv(zu_2 > z_bot1+dz,:)];
      z12 = transpose(z_top1:dz:zu_2(end));
      % M2 line 2
      [~,asv] = uvrot(u_2,v_2,phi_23);
      asv22 = [cumsum([dv2; asv(zu_2 == z_bot2+dz,:)],'reverse'); ...
        asv(zu_2 > z_bot2+dz,:)];
      z22 = transpose(z_top2:dz:zu_2(end));
      % M3
      [~,asv] = uvrot(u_3,v_3,phi_23);
      asv23 = [cumsum([dv2; asv(zu_3 == z_bot2+dz,:)],'reverse'); ...
        asv(zu_3 > z_bot2+dz,:)];
      z23 = transpose(z_top2:dz:zu_3(end));
      
      
      %% Apply different surface and cross-channel interpolations
      % Initialize ASVs
      ASV = ASV_G;
      ASV11 = squeeze(ASV_G(:,i_d1,:));
      ASV12 = squeeze(ASV_G(:,i_d2-1,:));
      ASV22 = squeeze(ASV_G(:,i_d2,:));
      ASV23 = squeeze(ASV_G(:,i_d3,:));
      
      % Insert geostrophic asv into full depth matix
      ASV11(z >= z11(1) & z <= z11(end),:) = asv11;
      ASV12(z >= z12(1) & z <= z12(end),:) = asv12;
      ASV22(z >= z22(1) & z <= z22(end),:) = asv22;
      ASV23(z >= z23(1) & z <= z23(end),:) = asv23;
      
      % For bottom referenced vel
      % Get index of first velocity below the surface
      i11 = find(z == z11(1));    i12 = find(z == z12(1));
      i22 = find(z == z22(1));    i23 = find(z == z23(1));
      switch i_surf
        case 1
          % Apply linear criteria
          ASV11(1:i11+1,:) = interp1(z(i11:i11+1),...
            ASV11(i11:i11+1,:),z(1:i11+1),'linear','extrap');
          ASV12(1:i12+1,:) = interp1(z(i12:i12+1),...
            ASV12(i12:i12+1,:),z(1:i12+1),'linear','extrap');
          ASV22(1:i22+1,:) = interp1(z(i22:i22+1),...
            ASV22(i22:i22+1,:),z(1:i22+1),'linear','extrap');
          ASV23(1:i23+1,:) = interp1(z(i23:i23+1),...
            ASV23(i23:i23+1,:),z(1:i23+1),'linear','extrap');
        case 2
          % Apply slab criteria
          ASV11(1:i11-1,:) = repmat(ASV11(i11,:),i11-1,1);
          ASV12(1:i12-1,:) = repmat(ASV12(i12,:),i12-1,1);
          ASV22(1:i22-1,:) = repmat(ASV22(i22,:),i22-1,1);
          ASV23(1:i23-1,:) = repmat(ASV23(i23,:),i23-1,1);
        case 3
          % Linearly interpolate to surface geostrophic
          % Insert surface geo to top
          [~,ASV11(1,:)] = uvrot(ug_1(1,:),vg_1(1,:),phi_12);
          [~,ASV12(1,:)] = uvrot(ug_2(1,:),vg_2(1,:),phi_12);
          [~,ASV22(1,:)] = uvrot(ug_2(1,:),vg_2(1,:),phi_23);
          [~,ASV23(1,:)] = uvrot(ug_3(1,:),vg_3(1,:),phi_23);
          % Interpolate
          ASV11(1:i11,:) = interp1([z(1); z(i11)],...
            [ASV11(1,:); ASV11(i11,:)],z(1:i11));
          ASV12(1:i12,:) = interp1([z(1); z(i12)],...
            [ASV12(1,:); ASV12(i12,:)],z(1:i12));
          ASV22(1:i22,:) = interp1([z(1); z(i22)],...
            [ASV22(1,:); ASV22(i22,:)],z(1:i22));
          ASV23(1:i23,:) = interp1([z(1); z(i23)],...
            [ASV23(1,:); ASV23(i23,:)],z(1:i23));
          % Add in Ekman component
          [~,asvE] = uvrot(ue_1,ve_1,phi_12);
          ASV11(z < 1000,:) = ASV11(z < 1000,:) + ...
            asvE(zu_1 < 1000,:);
          [~,asvE] = uvrot(ue_2,ve_2,phi_12);
          ASV12(z < 1000,:) = ASV12(z < 1000,:) + ...
            asvE(zu_2 < 1000,:);
          [~,asvE] = uvrot(ue_2,ve_2,phi_23);
          ASV22(z < 1000,:) = ASV22(z < 1000,:) + ...
            asvE(zu_2 < 1000,:);
          [~,asvE] = uvrot(ue_3,ve_3,phi_23);
          ASV23(z < 1000,:) = ASV23(z < 1000,:) + ...
            asvE(zu_3 < 1000,:);
        case 4
          % Linearly interpolate to surface currents
          % Insert surface vel to top
          [~,ASV11(1,:)] = uvrot(us_1,vs_1,phi_12);
          [~,ASV12(1,:)] = uvrot(us_2,vs_2,phi_12);
          [~,ASV22(1,:)] = uvrot(us_2,vs_2,phi_23);
          [~,ASV23(1,:)] = uvrot(us_3,vs_3,phi_23);
          % Interpolate
          ASV11(1:i11,:) = interp1([z(1); z(i11)],...
            [ASV11(1,:); ASV11(i11,:)],z(1:i11));
          ASV12(1:i12,:) = interp1([z(1); z(i12)],...
            [ASV12(1,:); ASV12(i12,:)],z(1:i12));
          ASV22(1:i22,:) = interp1([z(1); z(i22)],...
            [ASV22(1,:); ASV22(i22,:)],z(1:i22));
          ASV23(1:i23,:) = interp1([z(1); z(i23)],...
            [ASV23(1,:); ASV23(i23,:)],z(1:i23));
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
      
      % Insert asv profiles into full channel grid
      ASV(:,i_d1,:) = ASV11;
      ASV(:,i_d2-1,:) = ASV12;
      ASV(:,i_d2,:) = ASV22;
      ASV(:,i_d3,:) = ASV23;
      
      %% Interpolate between moorings and extrapolate to side walls
      [ASV,~] = channelInterpSS(ASV,ASV,bathy.dist,...
        [i_d1 i_d2-1 i_d2 i_d3],xpass);
      
      
      %% Interpolate z_sgth27 across channel
      z_sgth267 = [z1_sgth267' z1_sgth267' z2_sgth267' z3_sgth267' ...
        z3_sgth267'];
      z_sgth267 = interp1(bathy.dist([1 i_d1 i_d2 i_d3 end]),...
        z_sgth267',bathy.dist)';
      z_sgth275 = [z1_sgth275' z1_sgth275' z2_sgth275' z3_sgth275' ...
        z3_sgth275'];
      z_sgth275 = interp1(bathy.dist([1 i_d1 i_d2 i_d3 end]),...
        z_sgth275',bathy.dist)';
      
      
      %% Plot data
      
      % save xpassage scheme
      params_1.xpassage_method = ['geo' Gref 'Ref_' xpass];
      params_1.bottom_method = bot;
      params_1.surface_method = surf;
      
      % Collect instrument mean locations
      inst_xz = [bathy.dist(i_d1)*...
        ones(size(nanmean(RHO_1.intermediate.pres,2))) ...
        nanmean(RHO_1.intermediate.pres,2);...
        bathy.dist(i_d2)*...
        ones(size(nanmean(RHO_21.intermediate.pres,2))) ...
        nanmean(RHO_21.intermediate.pres,2);...
%         bathy.dist(i_d2)*...
%         ones(size(nanmean(RHO_22.intermediate.pres,2))) ...
%         nanmean(RHO_22.intermediate.pres,2);...
        bathy.dist(i_d3)*...
        ones(size(nanmean(RHO_3.intermediate.pres,2))) ...
        nanmean(RHO_3.intermediate.pres,2)];
      
      % Plot
      channelPlot_G(ASV,z,time,bathy,params_1,inst_xz,RHO,TEM,SAL)
      close all
      
      
      %% Save data
      i_moor = [i_d1 i_d2-1 i_d2 i_d3];
      f_out = [driveName dataPath params_1.channel '_' ...
        params_1.xpassage_method 'Xsection_' bot '_' surf '.mat'];
      save(f_out,'ASV','z','time','bathy','z_sgth267','z_sgth275',...
        'RHO','TEM','SAL','inst_xz','i_moor')
      
      
    end
  end
end