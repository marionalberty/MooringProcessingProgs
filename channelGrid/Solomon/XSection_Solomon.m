% XSection_Solomon.m

clear all; close all; clc

set(0,'defaultaxesfontsize',12,'defaultaxeslinewidth',0.7,...
  'defaultlinelinewidth',1,'defaultpatchlinewidth',0.7,...
  'defaultFigureColor','white')
%% Set path and directories
% Personnal paths:
driveName='/Users/marionsofiaalberty/MATLAB/Solomon_Sea/';

% Inputs/Outputs directories:
pathin =[driveName 'Moorings/Data/Gridded/SolomonStrait/'];
figout =[driveName 'Moorings/Figures/Gridded/SolomonStrait/'];


%% Load mooring lat & lon
% M1
load([pathin 'Solomon_M1/Solomon_M1.mat'],'params')
mlat(1) = params.lat;
mlon(1) = params.lon;
mbot(1) = params.bottomDepth;

% M2b
load([pathin 'Solomon_M2b/Solomon_M2b.mat'],'params')
mlat(2) = params.lat;
mlon(2) = params.lon;
mbot(2) = params.bottomDepth;

% M3
load([pathin 'Solomon_M3/Solomon_M3.mat'],'params')
mlat(3) = params.lat;
mlon(3) = params.lon;
mbot(3) = params.bottomDepth;

clear params


%% Get lat lon line for cross-section
p1 = polyfit(mlon(1:2),mlat(1:2),1);
p2 = polyfit(mlon(2:3),mlat(2:3),1);
% p3 = polyfit([mlon(1) mlon(3)],[mlat(1) mlat(3)],1);
p3 = polyfit(mlon,mlat,1);

% Pick out range of line
x1 = [152.903 mlon(2)];
y1 = polyval(p1,x1);
n1 = round(m_lldist(x1,y1)/2);

x2 = [mlon(2) 154.555];
y2 = polyval(p2,x2);
n2 = round(m_lldist(x2,y2)/2);

x3 = [152 154.555];
y3 = polyval(p3,x3);


%% Load different bathy sets
% GEBCO 2104 30 sec
[bathy_GEBCO.lon,bathy_GEBCO.lat] = meshgrid(...
  ncread('Solomon_Sea/GEBCO_2014/GEBCO_2014_SolomonSea.nc','lon')',...
  ncread('Solomon_Sea/GEBCO_2014/GEBCO_2014_SolomonSea.nc','lat')');
bathy_GEBCO.z = ncread('Solomon_Sea/GEBCO_2014/GEBCO_2014_SolomonSea.nc'...
  ,'elevation')';
bathy_GEBCO.SID = ...
  ncread('Solomon_Sea/GEBCO_2014/GEBCO_2014_SID_SolomonSea.nc','sid')';

% Sandwell 30 sec
load('Solomon_Sea/Bathy/topo30_sandwell.mat')
[bathy_Sandwell.lon,bathy_Sandwell.lat] = meshgrid(lon',lat');
bathy_Sandwell.z = double(z);
load('Solomon_Sea/Bathy/topo30_SID_sandwell.mat','SID')
bathy_Sandwell.SID = double(SID);
clear lat lon z SID


%% Extract bathy cross-section for each dataset
% Query points
xq1 = linspace(x1(1),x1(2),n1);
yq1 = polyval(p1,xq1);

xq2 = linspace(x2(1),x2(2),n2);
yq2 = polyval(p2,xq2);

xq = [xq1 xq2];     yq = [yq1 yq2];

% GEBCO
bathy_GEBCO.z_xsec = interp2(bathy_GEBCO.lon,bathy_GEBCO.lat,...
  bathy_GEBCO.z,xq,yq);
bathy_GEBCO.SID_xsec = interp2(bathy_GEBCO.lon,bathy_GEBCO.lat,...
  bathy_GEBCO.SID,xq,yq,'nearest');
bathy_GEBCO.lon_xsec = xq;
bathy_GEBCO.lat_xsec = yq;
bathy_GEBCO.dist_xsec = [0 cumsum(m_lldist(bathy_GEBCO.lon_xsec,...
  bathy_GEBCO.lat_xsec)')];

% Sandwell
bathy_Sandwell.z_xsec = interp2(bathy_Sandwell.lon,bathy_Sandwell.lat,...
  bathy_Sandwell.z,xq,yq);
bathy_Sandwell.SID_xsec = interp2(bathy_Sandwell.lon,bathy_Sandwell.lat,...
  bathy_Sandwell.SID,xq,yq,'nearest');
bathy_Sandwell.lon_xsec = xq;
bathy_Sandwell.lat_xsec = yq;
bathy_Sandwell.dist_xsec = [0 cumsum(m_lldist(bathy_Sandwell.lon_xsec,...
  bathy_Sandwell.lat_xsec)')];

%% Plot cross-section
figure
subplot(3,1,1:2)
plot(bathy_GEBCO.dist_xsec,bathy_GEBCO.z_xsec,'r')
hold on
plot(bathy_Sandwell.dist_xsec,bathy_Sandwell.z_xsec,'b')
x_cgrid = 0:10:200;
plot([x_cgrid; x_cgrid],[-4500*ones(size(x_cgrid)); ...
  zeros(size(x_cgrid))],'k')
x_mgrid = cumsum(m_lldist([x1(1) mlon],[y1(1) mlat]));
scatter(x_mgrid,zeros(size(x_mgrid)),100,'vr','fill')
scatter(x_mgrid,-mbot,100,'^r','fill')
legend('GEBCO','Sandwell','Channel Gridding','location',...
  'southwest')
axis tight
ylabel('Depth [m]')
xlabel('Distance across channel [km]')
title('Solomon Strait')
subplot(3,1,3)
plot(bathy_GEBCO.dist_xsec,bathy_GEBCO.SID_xsec ~= 0,'r')
hold on
plot(bathy_Sandwell.dist_xsec,bathy_Sandwell.SID_xsec ~= 0,'b')
axis tight
set(gca,'YLim',[0 1],'YTick',0:1,'YTickLabel',{'none','sounding'})

print([figout 'Xsection_bathyCompare.png'],'-dpng')


%% Load SADCP for velocity cross-section
sadcp.z = [10:20:990]';
sadcp.u = [];
sadcp.v = [];
sadcp.time = [];
sadcp.lat = [];
sadcp.lon = [];
% Pandora
load('Solomon_Sea/Pandora2012/SADCP/sadcp_38.mat')
sadcp.u = interp1(z_sadcp,u_sadcp,sadcp.z);
sadcp.v = interp1(z_sadcp,v_sadcp,sadcp.z);
[yy,mm,dd,HH,MM,SS] = jd2date(tim_sadcp);
sadcp.time = datenum(yy,mm,dd,HH,MM,SS)';
sadcp.lat = lat_sadcp';
sadcp.lon = lon_sadcp';

clear lat_sadcp lon_sadcp u_sadcp v_sadcp tim_sadcp z_sadcp

% MoorSPICE
load('Solomon_Sea/moorSPICE2014/SADCP/sadcp75nb.mat')
sadcp.u = [sadcp.u interp1(z(:,1),u,sadcp.z)];
sadcp.v = [sadcp.v interp1(z(:,1),v,sadcp.z)];
sadcp.time = [sadcp.time datenum'];
sadcp.lat = [sadcp.lat lat'];
sadcp.lon = [sadcp.lon lon'];

clear lat lon u v z datenum

% Archive data
path_sadcp = 'Solomon_Sea/sadcp_archive/';
arch = load([path_sadcp 'region_166_145_m15_m3_adcp_data.mat'],'sadcp');
for i = 1:numel(arch.sadcp)
  sadcp.u = [sadcp.u interp1(arch.sadcp(i).z,arch.sadcp(i).u./100,...
    sadcp.z)];
  sadcp.v = [sadcp.v interp1(arch.sadcp(i).z,arch.sadcp(i).v./100,...
    sadcp.z)];
  sadcp.time = [sadcp.time arch.sadcp(i).date(:)'];
  sadcp.lat = [sadcp.lat arch.sadcp(i).lat(:)'];
  sadcp.lon = [sadcp.lon arch.sadcp(i).lon(:)'];
end
% Remove missing data
null = 1.0000e+38;
sadcp.lat(sadcp.lon == null) = [];
sadcp.time(sadcp.lon == null) = [];
sadcp.u(:,sadcp.lon == null) = [];
sadcp.v(:,sadcp.lon == null) = [];
sadcp.lon(sadcp.lon == null) = [];


%% Indentify all data that isn't near the cross-section
% Max buffer between SADCP and line
d_buff = 20;
xq_v1 = linspace(x1(1),x1(2),numel(x_cgrid)/3);
yq_v1 = polyval(p1,xq_v1);
xq_v2 = linspace(x2(1),x2(2),2*numel(x_cgrid)/3);
yq_v2 = polyval(p2,xq_v2);
xq_v = [xq_v1 xq_v2(2:end)];
yq_v = [yq_v1 yq_v2(2:end)];

i_keep = nan(size(sadcp.lat));
i_q = i_keep;

for i = 1:numel(i_keep)
  dist = nan(size(xq_v));
  for k = 1:numel(dist)
    dist(k) = m_lldist([sadcp.lon(i) xq_v(k)],[sadcp.lat(i) yq_v(k)]); 
  end
  kk = sum(dist <= d_buff);
  if kk > 0
    % SADCP is within buffer
    i_keep(i) = 1;
    [~,i_q(i)] = min(dist);
  else
    % SADCP is not within buffer
    i_keep(i) = 0;
  end
end

clear dist

%% Keep data near line
sadcp.u = sadcp.u(:,logical(i_keep));
sadcp.v = sadcp.v(:,logical(i_keep));
sadcp.time = sadcp.time(logical(i_keep));
sadcp.lat = sadcp.lat(logical(i_keep));
sadcp.lon = sadcp.lon(logical(i_keep));
i_q = i_q(logical(i_keep));

clear i_keep kk i k

%% Average U & V along line
% get index of when the line switches
i_sw = numel(xq_v1);
% prelocate u and v
uq = nan(numel(sadcp.z),numel(xq_v));
vq = uq;
N_av = uq;

for k = 1:numel(xq_v)
  uq(:,k) = nanmean(sadcp.u(:,i_q == k),2);
  vq(:,k) = nanmean(sadcp.v(:,i_q == k),2);
  N_av(:,k) = sum(~isnan(sadcp.u(:,i_q == k)),2);
end
% Mask averages with fewer than 4 points
Nmask = ones(size(uq));
Nmask(N_av < 4) = nan;
uq = uq.*Nmask;
vq = vq.*Nmask;
% Prelocate rotated velocities matrix
acv = nan(numel(sadcp.z),numel(xq_v)+1);
xcv = acv;
% Rotate for first line
[xcv(:,1:i_sw),acv(:,1:i_sw)] = uvrot(uq(:,1:i_sw),vq(:,1:i_sw),...
  90-azimuth(mlat(1),mlon(1),mlat(2),mlon(2)));
% Rotate for second line
[xcv(:,i_sw+1:end),acv(:,i_sw+1:end)] = uvrot(uq(:,i_sw:end),...
  vq(:,i_sw:end),90-azimuth(mlat(2),mlon(2),mlat(3),mlon(3)));
% Calc distance along x-section
dist = [0 cumsum(m_lldist(xq_v([1:i_sw i_sw:end]),...
  yq_v([1:i_sw i_sw:end])))'];


%% Save SADCP mean section
data.xdist = dist;
data.lat = yq_v([1:i_sw i_sw:end]);
data.lon = xq_v([1:i_sw i_sw:end]);
data.z = sadcp.z;
data.asv = acv;
data.xsv = xcv;

save([driveName 'sadcp_archive/SolomonStrait_meanAsvXsv.mat'],'data')


%% Plot cross-sections of velocity
figure('position',[0 0 600 800])
subplot(211)
pcolor(dist,sadcp.z,acv)
shading interp
hold on
scatter(x_mgrid,zeros(size(x_mgrid)),100,'vk','fill')
scatter(x_mgrid,mbot,100,'^k','fill')
plot(bathy_Sandwell.dist_xsec,-bathy_Sandwell.z_xsec,'k')
axis ij
ylim([0 1000])
xlim([0 bathy_Sandwell.dist_xsec(end)])
colorbar
colormap(redblue(25))
caxis([-1 1]*max(abs(acv(:))))
xlabel('Distance Along Cross-Section [km]')
ylabel('Depth [m]')
title('Solomon Along Channel Velocity [m/s]')

subplot(212)
pcolor(dist,sadcp.z,xcv)
shading interp
hold on
scatter(x_mgrid,zeros(size(x_mgrid)),100,'vk','fill')
scatter(x_mgrid,mbot,100,'^k','fill')
plot(bathy_Sandwell.dist_xsec,-bathy_Sandwell.z_xsec,'k')
axis ij
ylim([0 1000])
xlim([0 bathy_Sandwell.dist_xsec(end)])
colorbar
colormap(redblue(25))
caxis([-1 1]*max(abs(xcv(:))))
xlabel('Distance Along Cross-Section [km]')
ylabel('Depth [m]')
title('Solomon Across Channel Velocity [m/s]')

print([figout 'Xsection_SADCPvelocity.png'],'-dpng')


%% Plot histogram of SADCP monthly distribution
[~,MM,~,~,~,~] = datevec(sadcp.time);
figure
histogram(MM,0.5:12.5,'normalization','probability')
xlabel('Month')
ylabel('% of Total Profiles')
title('Solomon Strait SADCP Monthly Breakdown')

print([figout 'Xsection_MonthHistogram.png'],'-dpng')


%% Calculate vertical shear from mean u & v
% Get dz
dz = mean(diff(sadcp.z));
% Get new z coords
z_shear = sadcp.z(1:end-1)+dz/2;
% calc vertical shear
davdz = diff(acv,1,1)./dz;
dxvdz = diff(xcv,1,1)./dz;


%% Plot vertical Shear
figure('position',[0 0 600 800])
subplot(211)
pcolor(dist,z_shear,davdz)
shading interp
hold on
scatter(x_mgrid,zeros(size(x_mgrid)),100,'vk','fill')
scatter(x_mgrid,mbot,100,'^k','fill')
plot(bathy_Sandwell.dist_xsec,-bathy_Sandwell.z_xsec,'k')
axis ij
ylim([0 1000])
xlim([0 bathy_Sandwell.dist_xsec(end)])
colorbar
colormap(redblue(25))
caxis([-1 1]*max(abs(davdz(:)))*0.5)
xlabel('Distance Along Cross-Section [km]')
ylabel('Depth [m]')
title('Solomon Along Channel Vertical Shear [1/s]')

subplot(212)
pcolor(dist,z_shear,dxvdz)
shading interp
hold on
scatter(x_mgrid,zeros(size(x_mgrid)),100,'vk','fill')
scatter(x_mgrid,mbot,100,'^k','fill')
plot(bathy_Sandwell.dist_xsec,-bathy_Sandwell.z_xsec,'k')
axis ij
ylim([0 1000])
xlim([0 bathy_Sandwell.dist_xsec(end)])
colorbar
colormap(redblue(25))
caxis([-1 1]*max(abs(dxvdz(:)))*0.5)
xlabel('Distance Along Cross-Section [km]')
ylabel('Depth [m]')
title('Solomon Across Channel Vertical Shear [1/s]')

print([figout 'Xsection_SADCPverticalShear.png'],'-dpng')


%% Plot mean profile of vertical shear
davdz_bar = mean(100*davdz,2,'omitnan');
davdz_std = std(100*davdz,0,2,'omitnan');
figure('position',[0 0 600 400])
herrorbar(davdz_bar,z_shear,davdz_std,'k.')
hold on
plot(davdz_bar,z_shear,'k')
axis ij
ylim([0 1000])
xlim([min(davdz_bar-davdz_std) max(davdz_bar+davdz_std)])
ylabel('Depth [m]')
title('Solomon Mean Across Channel Vertical Shear [cm s^{-1} m^{-1}]')

print([figout 'Xsection_SADCPverticalShearBar.png'],'-dpng')


%% Calculate horizontal shear from mean u & v
% get dx
dist_shear = dist(1:end-1) + diff(dist)./2;
dx = mean(diff(dist));
% calc shear
davdx = diff(acv,1,2)./dx;
dxvdx = diff(xcv,1,2)./dx;


%% Plot Horizontal Shear
figure('position',[0 0 600 800])
subplot(211)
pcolor(dist_shear,sadcp.z,davdx)
shading interp
hold on
scatter(x_mgrid,zeros(size(x_mgrid)),100,'vk','fill')
scatter(x_mgrid,mbot,100,'^k','fill')
plot(bathy_Sandwell.dist_xsec,-bathy_Sandwell.z_xsec,'k')
axis ij
ylim([0 1000])
xlim([0 bathy_Sandwell.dist_xsec(end)])
colorbar
colormap(redblue(25))
caxis([-1 1]*max(abs(davdx(:))))
xlabel('Distance Along Cross-Section [km]')
ylabel('Depth [m]')
title('Solomon Along Channel Horizontal Shear [m/s km]')

subplot(212)
pcolor(dist_shear,sadcp.z,dxvdx)
shading interp
hold on
scatter(x_mgrid,zeros(size(x_mgrid)),100,'vk','fill')
scatter(x_mgrid,mbot,100,'^k','fill')
plot(bathy_Sandwell.dist_xsec,-bathy_Sandwell.z_xsec,'k')
axis ij
ylim([0 1000])
xlim([0 bathy_Sandwell.dist_xsec(end)])
colorbar
colormap(redblue(25))
caxis([-1 1]*max(abs(dxvdx(:))))
xlabel('Distance Along Cross-Section [km]')
ylabel('Depth [m]')
title('Solomon Across Channel Horizontal Shear [m/s km]')

print([figout 'Xsection_SADCPhorizontalShear.png'],'-dpng')


%% Plot locations of grid lines compared to the mooring line
figure('position',[10 10 800 300])
m_proj('lambert','long',[152.8 154.6],'lat',[-5.3 -4.7])
m_contour(bathy_Sandwell.lon,bathy_Sandwell.lat,bathy_Sandwell.z,-4500:250:0)
hold on
m_contour(bathy_Sandwell.lon,bathy_Sandwell.lat,bathy_Sandwell.z,[0 0],'k',...
  'LineWidth',2)
m_plot([x1 x2],[y1 y2],'r')
m_scatter(sadcp.lon,sadcp.lat,10,'k','fill')
m_scatter(mlon,mlat,300,'rp','fill')
m_grid('box','fancy','tickdir','in');

print([figout 'Xsection_map.png'],'-dpng')


%% Save bathy cross-section
% Query points
xq1 = linspace(x1(1),x1(2),n1);
yq1 = polyval(p1,xq1);

xq2 = linspace(x2(1),x2(2),n2);
yq2 = polyval(p2,xq2);

xq = [xq1 xq2];     yq = [yq1 yq2];

% Sandwell
bathy.z = -1*interp2(bathy_Sandwell.lon,bathy_Sandwell.lat,...
  bathy_Sandwell.z,xq,yq);
bathy.SID = interp2(bathy_Sandwell.lon,bathy_Sandwell.lat,...
  bathy_Sandwell.SID,xq,yq,'nearest');
bathy.lon = xq;
bathy.lat = yq;
bathy.dist = [0 cumsum(m_lldist(xq,yq)')];

% save
save([pathin 'Xsection_Bathy.mat'],'bathy')
