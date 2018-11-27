% mk_SStimeseriesG.m

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
g = -9.8;   %[m/s^2]


%% Load bathy and initialize channel grid
% bathy
load([driveName dataPath 'Xsection_Bathy.mat'])

% z-grid
z = transpose(10:dz:max(bathy.z));


%% Load mooring data
% Solomon M1
load([driveName dataPath 'Solomon_M1/Solomon_M1.mat'])
UV_1 = UV;    SS_1 = SS;      TT_1 = TT;    params_1 = params;
load([driveName dataPath ...
  'Solomon_M1/Solomon_M1_slab2bot_sat4surf.mat'],'UV')
UV_1.U_geo = UV.U_geo;          UV_1.V_geo = UV.V_geo;
UV_1.U_ekm = UV.U - UV.U_geo;   UV_1.V_ekm = UV.V - UV.V_geo;

% Solomon M2a and Solomon M2b
load([driveName dataPath 'Solomon_M2a/Solomon_M2a.mat'])
SS_2 = SS;      TT_2 = TT;
load([driveName dataPath 'Solomon_M2b/Solomon_M2b.mat'])
UV_2 = UV;    params_2 = params;
load([driveName dataPath ...
  'Solomon_M2b/Solomon_M2b_slab2bot_sat4surf.mat'],'UV')
UV_2.U_geo = UV.U_geo;          UV_2.V_geo = UV.V_geo;
UV_2.U_ekm = UV.U - UV.U_geo;   UV_2.V_ekm = UV.V - UV.V_geo;

% Solomon M3
load([driveName dataPath 'Solomon_M3/Solomon_M3.mat'])
UV_3 = UV;    SS_3 = SS;      TT_3 = TT;    params_3 = params;
load([driveName dataPath ...
  'Solomon_M3/Solomon_M3_slab2bot_sat4surf.mat'],'UV')
UV_3.U_geo = UV.U_geo;          UV_3.V_geo = UV.V_geo;
UV_3.U_ekm = UV.U - UV.U_geo;   UV_3.V_ekm = UV.V - UV.V_geo;

% Housekeeping
clear UV TT SS params


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


%% Filer and subsample density data
% Solomon M1
% Extract data to be filtered
M1.temp = naninterp(TT_1.temp')';   M1.psal = naninterp(SS_1.psal')';
M1.u = naninterp(UV_1.U')';         M1.v = naninterp(UV_1.V')';
M1.ug = naninterp(UV_1.U_geo')';    M1.vg = naninterp(UV_1.V_geo')';
M1.ue = naninterp(UV_1.U_ekm')';    M1.ve = naninterp(UV_1.V_ekm')';
% Filter
M_filt = moorFilter(M1,TT_1.time,tfilt);
% Extract full coverage profiles
i_s = find(sum(isnan(M_filt.psal),2) == 0);
i_u = find(sum(isnan(UV_1.U),2)/numel(UV_1.time) < 0.25);
% Subsample
t_1 = interp1(TT_1.time,M_filt.temp(i_s,:)',time)';
s_1 = interp1(TT_1.time,M_filt.psal(i_s,:)',time)';
z_1 = TT_1.z(i_s);
u_1 = interp1(UV_1.time,M_filt.u(i_u,:)',time)';
v_1 = interp1(UV_1.time,M_filt.v(i_u,:)',time)';
ug_1 = interp1(UV_1.time,M_filt.ug',time)';
vg_1 = interp1(UV_1.time,M_filt.vg',time)';
ue_1 = interp1(UV_1.time,M_filt.ue',time)';
ve_1 = interp1(UV_1.time,M_filt.ve',time)';
zu_1 = UV_1.z(i_u);
zug_1 = UV_1.z;

% Solomon M2
% Extract data to be filtered
M2.temp = naninterp(TT_2.temp')';   M2.psal = naninterp(SS_2.psal')';
M2.u = naninterp(UV_2.U')';         M2.v = naninterp(UV_2.V')';
M2.ug = naninterp(UV_2.U_geo')';    M2.vg = naninterp(UV_2.V_geo')';
M2.ue = naninterp(UV_2.U_ekm')';    M2.ve = naninterp(UV_2.V_ekm')';
% Filter
M_filt = moorFilter(M2,TT_2.time,tfilt);
% Extract full coverage profiles
i_s = find(sum(isnan(M_filt.psal),2) == 0);
i_u = find(sum(isnan(UV_2.U),2)/numel(UV_2.time) < 0.25);
% Subsample
t_2 = interp1(TT_2.time,M_filt.temp(i_s,:)',time)';
s_2 = interp1(TT_2.time,M_filt.psal(i_s,:)',time)';
z_2 = TT_2.z(i_s);
u_2 = interp1(UV_2.time,M_filt.u(i_u,:)',time)';
v_2 = interp1(UV_2.time,M_filt.v(i_u,:)',time)';
ug_2 = interp1(UV_2.time,M_filt.ug',time)';
vg_2 = interp1(UV_2.time,M_filt.vg',time)';
ue_2 = interp1(UV_2.time,M_filt.ue',time)';
ve_2 = interp1(UV_2.time,M_filt.ve',time)';
zu_2 = UV_2.z(i_u);
zug_2 = UV_2.z;

% Solomon M3
% Extract data to be filtered
M3.temp = naninterp(TT_3.temp')';   M3.psal = naninterp(SS_3.psal')';
M3.u = naninterp(UV_3.U')';         M3.v = naninterp(UV_3.V')';
M3.ug = naninterp(UV_3.U_geo')';    M3.vg = naninterp(UV_3.V_geo')';
M3.ue = naninterp(UV_3.U_ekm')';    M3.ve = naninterp(UV_3.V_ekm')';
% Filter
M_filt = moorFilter(M3,TT_3.time,tfilt);
% Extract full coverage profiles
i_s = find(sum(isnan(M_filt.psal),2) == 0);
i_u = find(sum(isnan(UV_3.U),2)/numel(UV_3.time) < 0.25);
% Subsample
t_3 = interp1(TT_3.time,M_filt.temp(i_s,:)',time)';
s_3 = interp1(TT_3.time,M_filt.psal(i_s,:)',time)';
z_3 = TT_3.z(i_s);
u_3 = interp1(UV_3.time,M_filt.u(i_u,:)',time)';
v_3 = interp1(UV_3.time,M_filt.v(i_u,:)',time)';
ug_3 = interp1(UV_3.time,M_filt.ug',time)';
vg_3 = interp1(UV_3.time,M_filt.vg',time)';
ue_3 = interp1(UV_3.time,M_filt.ue',time)';
ve_3 = interp1(UV_3.time,M_filt.ve',time)';
zu_3 = UV_3.z(i_u);
zug_3 = UV_3.z;

% Housekeeping
clear M* i_u i_s TT_* UV_*


%% Calculate density at each mooring
rho_1 = sw_dens(s_1,t_1,0);
rho_2 = sw_dens(s_2,t_2,0);
rho_3 = sw_dens(s_3,t_3,0);


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


%% Calculate vertical profile of shear for each line
for i_Gref = 1:2
  % West Segment (M1 & M2, line 1)
  % drho/dx
  z_top1 = max([z_1(1) z_2(1)]);
  z_bot1 = min([z_1(end) z_2(end) bathy.z(i_d1:i_d2)]);
  i_zt = find(z == z_top1);
  i_zb = find(z <= z_bot1,1,'last');
  z_bot1 = z(i_zb);
  drho1 = rho_2(z_2 >= z_top1 & z_2 <= z_bot1,:) - ...
    rho_1(z_1 >= z_top1 & z_1 <= z_bot1,:);
  dx1 = (bathy.dist(i_d2) - bathy.dist(i_d1))*1000;
  % grf = g/rho_0 f
  grf1 = repmat(g/sw_f(mean(bathy.lat(i_d1:i_d2)))./mean(...
    [rho_2(z_2 >= z_top1 & z_2 <= z_bot1,:); ...
    rho_1(z_1 >= z_top1 & z_1 <= z_bot1,:)]),i_zb-i_zt+1,1);
  % dv
  dv1 = dz * grf1 .* drho1 ./ dx1;
  
  % East Segment (M2 & M3, line 2)
  % drho/dx
  z_top2 = max([z_2(1) z_3(1)]);
  z_bot2 = min([z_2(end) z_3(end)]);
  i_zt = find(z == z_top2);
  i_zb = find(z == z_bot2);
  drho2 = rho_3(z_3 >= z_top2 & z_3 <= z_bot2,:) - ...
    rho_2(z_2 >= z_top2 & z_2 <= z_bot2,:);
  dx2 = (bathy.dist(i_d3) - bathy.dist(i_d2))*1000;
  % grf = g/rho_0 f
  grf2 = repmat(g/sw_f(mean(bathy.lat(i_d2:i_d3)))./mean(...
    [rho_3(z_3 >= z_top2 & z_3 <= z_bot2,:); ...
    rho_2(z_2 >= z_top2 & z_2 <= z_bot2,:)]),i_zb-i_zt+1,1);
  % dv
  dv2 = dz * grf2 .* drho2 ./ dx2;
  
  switch i_Gref
    case 1
      % Use near surface velocities for reference
      % Set case for reference velocity
      Gref = 'Top';
      % M1
      [~,asv] = uvrot(ug_1,vg_1,phi_12);
      asv11 = [asv(zug_1 < z_top1-dz,:); ...
        cumsum([asv(zug_1 == z_top1-dz,:); -1*dv1])];
      z11 = transpose(zug_1(1):dz:z_bot1);
      % M2 line 1
      [~,asv] = uvrot(ug_2,vg_2,phi_12);
      asv12 = [asv(zug_2 < z_top1-dz,:); ...
        cumsum([asv(zug_2 == z_top1-dz,:); -1*dv1])];
      z12 = transpose(zug_2(1):dz:z_bot1);
      % Use near bottom velocities for reference
      % M2 line 2
      [~,asv] = uvrot(ug_2,vg_2,phi_23);
      asv22 = [asv(zug_2 < z_top2-dz,:); ...
        cumsum([asv(zug_2 == z_top2-dz,:); -1*dv2])];
      z22 = transpose(zug_2(1):dz:z_bot2);
      % M3
      [~,asv] = uvrot(ug_3,vg_3,phi_23);
      asv23 = [asv(zug_3 < z_top2-dz,:); ...
        cumsum([asv(zug_3 == z_top2-dz,:); -1*dv2])];
      z23 = transpose(zug_3(1):dz:z_bot2);
    case 2
      % Use near bottom velocities for reference
      % Set case for reference velocity
      Gref = 'Bot';
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
      % Use near bottom velocities for reference
      % M2 line 2
      [~,asv] = uvrot(u_2,v_2,phi_23);
      asv22 = [cumsum([dv2; asv(zu_2 == z_bot2+dz,:)],'reverse'); ...
        asv(zu_2 > z_bot2+dz,:)];
      z22 = transpose(z_top2:dz:zu_2(end));
      % M3
      [~,asv] = uvrot(u_3,v_3,phi_23);
      asv23 = cumsum([dv2; asv(end,:)],'reverse');
      z23 = transpose(z_top2:dz:z_bot2+dz);
  end
  
  
  %% Apply different surface, bottom and cross-channel interpolations
  for i_surf = 1:3
    for i_bot = 1:2
      for i_xpass = 1:3
        % Initialize ASVs
        ASV = ASV_G;
        ASV11 = squeeze(ASV_G(:,i_d1,:));
        ASV12 = squeeze(ASV_G(:,i_d2,:));
        ASV22 = squeeze(ASV_G(:,i_d2,:));
        ASV23 = squeeze(ASV_G(:,i_d3,:));
        
        % Insert geostrophic asv into full depth matix
        ASV11(z >= z11(1) & z <= z11(end),:) = asv11;
        ASV12(z >= z12(1) & z <= z12(end),:) = asv12;
        ASV22(z >= z22(1) & z <= z22(end),:) = asv22;
        ASV23(z >= z23(1) & z <= z23(end),:) = asv23;
        
        switch i_Gref
          %% Mooring surface extrapolation method
          case 1
            % For top referenced vel
            surf = 'sat4surf';
            % Add in Ekman component
            [~,asvE] = uvrot(ue_1,ve_1,phi_12);
            ASV11(z < zu_1(end),:) = ASV11(z < zu_1(end),:) + ...
              asvE(zug_1 < zu_1(end),:);
            [~,asvE] = uvrot(ue_2,ve_2,phi_12);
            ASV12(z < zu_2(end),:) = ASV12(z < zu_2(end),:) + ...
              asvE(zug_2 < zu_2(end),:);
            [~,asvE] = uvrot(ue_2,ve_2,phi_23);
            ASV22(z < zu_2(end),:) = ASV22(z < zu_2(end),:) + ...
              asvE(zug_2 < zu_2(end),:);
            [~,asvE] = uvrot(ue_3,ve_3,phi_23);
            ASV23(z < zu_3(end),:) = ASV23(z < zu_3(end),:) + ...
              asvE(zug_3 < zu_3(end),:);
            
          case 2
            % For bottom referenced vel
            % Get index of first velocity below the surface
            i11 = find(z == z11(1));    i12 = find(z == z12(1));
            i22 = find(z == z22(1));    i23 = find(z == z23(1));
            switch i_surf
              case 1
                surf = 'linear2surf';
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
                surf = 'slab2surf';
                % Apply slab criteria
                ASV11(1:i11-1,:) = repmat(ASV11(i11,:),i11-1,1);
                ASV12(1:i12-1,:) = repmat(ASV12(i12,:),i12-1,1);
                ASV22(1:i22-1,:) = repmat(ASV22(i22,:),i22-1,1);
                ASV23(1:i23-1,:) = repmat(ASV23(i23,:),i23-1,1);
              case 3
                surf = 'sat4surf';
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
                ASV11(z < zu_1(end),:) = ASV11(z < zu_1(end),:) + ...
                  asvE(zug_1 < zu_1(end),:);
                [~,asvE] = uvrot(ue_2,ve_2,phi_12);
                ASV12(z < zu_2(end),:) = ASV12(z < zu_2(end),:) + ...
                  asvE(zug_2 < zu_2(end),:);
                [~,asvE] = uvrot(ue_2,ve_2,phi_23);
                ASV22(z < zu_2(end),:) = ASV22(z < zu_2(end),:) + ...
                  asvE(zug_2 < zu_2(end),:);
                [~,asvE] = uvrot(ue_3,ve_3,phi_23);
                ASV23(z < zu_3(end),:) = ASV23(z < zu_3(end),:) + ...
                  asvE(zug_3 < zu_3(end),:);
            end
        end
        
        %% Mooring bottom extrapolation method
        % Get index of last velocity for interpolation
        i11 = find(z == z11(end));    i12 = find(z == z12(end));
        i22 = find(z == z22(end));    i23 = find(z == z23(end));
        % Get index of last bin to interpolate to before the bottom
        n1 = find(isnan(ASV11(:,1)),1)-2;
        n2 = find(isnan(ASV22(:,1)),1)-2;
        n3 = find(isnan(ASV23(:,1)),1)-2;
        switch i_bot
          case 1
            bot = 'linear2bot';
            % Apply linear criteria
            ASV11(i11:n1+1,:) = interp1(z([i11 n1+1]),...
              ASV11([i11 n1+1],:),z(i11:n1+1));
            ASV12(i12:n2+1,:) = interp1(z([i12 n2+1]),...
              ASV12([i12 n2+1],:),z(i12:n2+1));
            ASV22(i22:n2+1,:) = interp1(z([i22 n2+1]),...
              ASV22([i22 n2+1],:),z(i22:n2+1));
            ASV23(i23:n3+1,:) = interp1(z([i23 n3+1]),...
              ASV23([i23 n3+1],:),z(i23:n3+1));
          case 2
            bot = 'slab2bot';
            % Apply slab criteria
            ASV11(i11+1:n1,:) = repmat(ASV11(i11,:),n1-i11,1);
            ASV12(i12+1:n2,:) = repmat(ASV12(i12,:),n2-i12,1);
            ASV22(i22+1:n2,:) = repmat(ASV22(i22,:),n2-i22,1);
            ASV23(i23+1:n3,:) = repmat(ASV23(i23,:),n3-i23,1);
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
          ones(size(nanmean(SS_1.intermediate.pres,2))) ...
          nanmean(SS_1.intermediate.pres,2);...
          bathy.dist(i_d2)*...
          ones(size(nanmean(SS_2.intermediate.pres,2))) ...
          nanmean(SS_2.intermediate.pres,2);...
          bathy.dist(i_d3)*...
          ones(size(nanmean(SS_3.intermediate.pres,2))) ...
          nanmean(SS_3.intermediate.pres,2)];
        
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
end