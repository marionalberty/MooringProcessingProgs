% BathyXSection_Vitiaz.m

clear all; close all; clc

set(0,'defaultaxesfontsize',12,'defaultaxeslinewidth',0.7,...
  'defaultlinelinewidth',1,'defaultpatchlinewidth',0.7,...
  'defaultFigureColor','white')
%% Set path and directories

% Personnal paths:
driveName='/Users/marionsofiaalberty/MATLAB/Solomon_Sea/';

% Inputs/Outputs directories:
pathin =[driveName 'Moorings/Data/Gridded/Vitiaz/'];
figout =[driveName 'Moorings/Figures/Gridded/Vitiaz/'];


%% Load mooring lat & lon
% East
load([pathin 'VitiazEast/VitiazEast.mat'],'params')
mlat(3) = params.lat;
mlon(3) = params.lon;
mbot(3) = params.bottomDepth;

% Middle
load([pathin 'VitiazMiddle/VitiazMiddle.mat'],'params')
mlat(2) = params.lat;
mlon(2) = params.lon;
mbot(2) = params.bottomDepth;

% West
load([pathin 'VitiazWest/VitiazWest.mat'],'params')
mlat(1) = params.lat;
mlon(1) = params.lon;
mbot(1) = params.bottomDepth;

clear params


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


%% Get line for cross-section

% Lat Lon Limits to use while fitting line to coast
% [Lat1 Lon1 Lat2 Lon2]
% coast_lat = [ -5.6667  -5.7000  -6.3000  -6.1000];
% coast_lon = [147.6000 147.6000 148.3333 148.3333];
coast_lat = [ -6.0  -5.8  -6.3  -6.3];
coast_lon = [147.3 147.5 148.0 147.6];

% Get indecies of coast
% i_coast = find(bathy_GEBCO.z > -20 & bathy_GEBCO.z < 20 & ...
% bathy_GEBCO.lat < coast_lim(1) & bathy_GEBCO.lon > coast_lim(2) & ...
% bathy_GEBCO.lat > coast_lim(3) & bathy_GEBCO.lon < coast_lim(4));
% i_coast = find(bathy_GEBCO.z > -900 & bathy_GEBCO.z < -800 & ...
%   inpolygon(bathy_GEBCO.lon,bathy_GEBCO.lat,coast_lon,coast_lat));
i_coast = find(bathy_GEBCO.z > -250 & bathy_GEBCO.z < -150 & ...
  inpolygon(bathy_GEBCO.lon,bathy_GEBCO.lat,coast_lon,coast_lat));

% Get lat and lon for polyfit
lat_coast = bathy_GEBCO.lat(i_coast);
lon_coast = bathy_GEBCO.lon(i_coast);

% Get line that describes the coast
p_coast = polyfit(lon_coast,lat_coast,1);

% Use recipricol of coast slope for cross-section line
% Find line's y-intercept by running line through middle mooring
p(1) = -1/p_coast(1);
p(2) = mlat(2) - (p(1) * mlon(2));

% Pick out cross-section line coordinates
x = [147.62 147.94];
y = polyval(p,x);

% number of transport grid points
n_gp = round(m_lldist(x,y));


%% Extract bathy cross-section for each dataset

% Query point
xq = linspace(x(1),x(2),100);
yq = polyval(p,xq);

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

%% Plot cross-section of depth

% Get across-channel distance for each mooring
i_q = nan(size(mlon));
for i = 1:numel(mlon)
  dist = nan(size(xq));
  for k = 1:numel(dist)
    dist(k) = m_lldist([mlon(i) xq(k)],[mlat(i) yq(k)]);
  end
  [~,i_q(i)] = min(dist);
end
x_mgrid = bathy_GEBCO.dist_xsec(i_q);

% make figure
figure
subplot(3,1,1:2)
plot(bathy_GEBCO.dist_xsec,bathy_GEBCO.z_xsec,'r')
hold on
plot(bathy_Sandwell.dist_xsec,bathy_Sandwell.z_xsec,'b')
x_cgrid = linspace(0,m_lldist(x,y),n_gp);
plot([x_cgrid; x_cgrid],[-1200*ones(size(x_cgrid)); ...
  zeros(size(x_cgrid))],'k')
scatter(x_mgrid,zeros(size(x_mgrid)),100,'vr','fill')
scatter(x_mgrid,-mbot,100,'^r','fill')
legend('GEBCO','Sandwell','Channel Gridding','location',...
  'southeast')
axis tight
ylabel('Depth [m]')
xlabel('Distance across channel [km]')
title('Vitiaz Strait')
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
flist_sadcp = dirs(fullfile(path_sadcp,'*.nc'));
for i = 1:numel(flist_sadcp)
  sadcp.u = [sadcp.u interp1(ncread([path_sadcp flist_sadcp(i).name],...
    'depth'),ncread([path_sadcp flist_sadcp(i).name],...
    'u'),sadcp.z)];
  sadcp.v = [sadcp.v interp1(ncread([path_sadcp flist_sadcp(i).name],...
    'depth'),ncread([path_sadcp flist_sadcp(i).name],...
    'v'),sadcp.z)];
  sadcp.time = [sadcp.time ncread([path_sadcp flist_sadcp(i).name],...
    'decday')'+datenum(1992,1,1)];
  sadcp.lat = [sadcp.lat ncread([path_sadcp flist_sadcp(i).name],...
    'latitude')'];
  sadcp.lon = [sadcp.lon ncread([path_sadcp flist_sadcp(i).name],...
    'longitude')'];
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
xq_v = linspace(x(1),x(2),round(numel(x_cgrid)/2));
yq_v = polyval(p,xq_v);

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
% Rotate
[xcv,acv] = uvrot(uq,vq,90-azimuth(yq_v(1),xq_v(1),yq_v(end),xq_v(end)));
% Calc distance along x-section
dist = [0 cumsum(m_lldist(xq_v,yq_v))'];


%% Save SADCP mean section
data.xdist = dist;
data.lat = yq_v;
data.lon = xq_v;
data.z = sadcp.z;
data.asv = acv;
data.xsv = xcv;

save([driveName 'sadcp_archive/Vitiaz_meanAsvXsv.mat'],'data')


%% Plot cross-sections of velocity
figure('position',[0 0 600 800])
subplot(211)
pcolor(dist,sadcp.z,acv)
shading interp
hold on
scatter(x_mgrid,zeros(size(x_mgrid)),100,'vk','fill')
scatter(x_mgrid,mbot,100,'^k','fill')
plot(bathy_GEBCO.dist_xsec,-bathy_GEBCO.z_xsec,'k')
axis ij tight
colorbar
colormap(redblue(25))
caxis([-1 1]*max(abs(acv(:)))*0.75)
xlabel('Distance Along Cross-Section [km]')
ylabel('Depth [m]')
title('Vitiaz Along Channel Velocity [m/s]')

subplot(212)
pcolor(dist,sadcp.z,xcv)
shading interp
hold on
scatter(x_mgrid,zeros(size(x_mgrid)),100,'vk','fill')
scatter(x_mgrid,mbot,100,'^k','fill')
plot(bathy_GEBCO.dist_xsec,-bathy_GEBCO.z_xsec,'k')
axis ij tight
colorbar
colormap(redblue(25))
caxis([-1 1]*max(abs(xcv(:)))*0.75)
xlabel('Distance Along Cross-Section [km]')
ylabel('Depth [m]')
title('Vitiaz Across Channel Velocity [m/s]')

print([figout 'Xsection_SADCPvelocity.png'],'-dpng')


%% Plot histogram of SADCP monthly distribution
[~,MM,~,~,~,~] = datevec(sadcp.time);
figure
histogram(MM,0.5:12.5,'normalization','probability')
xlabel('Month')
ylabel('% of Total Profiles')
title('Vitiaz Strait SADCP Monthly Breakdown')

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
plot(bathy_GEBCO.dist_xsec,-bathy_GEBCO.z_xsec,'k')
axis ij tight
colorbar
colormap(redblue(25))
caxis([-1 1]*max(abs(davdz(:)))*0.3)
xlabel('Distance Along Cross-Section [km]')
ylabel('Depth [m]')
title('Vitiaz Along Channel Vertical Shear [1/s]')

subplot(212)
pcolor(dist,z_shear,dxvdz)
shading interp
hold on
scatter(x_mgrid,zeros(size(x_mgrid)),100,'vk','fill')
scatter(x_mgrid,mbot,100,'^k','fill')
plot(bathy_GEBCO.dist_xsec,-bathy_GEBCO.z_xsec,'k')
axis ij tight
colorbar
colormap(redblue(25))
caxis([-1 1]*max(abs(dxvdz(:))))
xlabel('Distance Along Cross-Section [km]')
ylabel('Depth [m]')
title('Vitiaz Across Channel Vertical Shear [1/s]')

print([figout 'Xsection_SADCPverticalShear.png'],'-dpng')


%% Plot mean profile of vertical shear
davdz_bar = mean(100*davdz,2,'omitnan');
davdz_std = std(100*davdz,0,2,'omitnan');
figure('position',[0 0 600 400])
herrorbar(davdz_bar,z_shear,davdz_std,'k.')
hold on
plot(davdz_bar,z_shear,'k')
axis ij
ylim([0 max(mbot)])
xlim([min(davdz_bar-davdz_std) max(davdz_bar+davdz_std)])
ylabel('Depth [m]')
title('Vitiaz Mean Across Channel Vertical Shear [cm s^{-1} m^{-1}]')

print([figout 'Xsection_SADCPverticalShearBar.png'],'-dpng')


%% Calculate horizontal shear from mean u & v
% get dx
dist_shear = dist(1:end-1) + diff(dist)./2;
dx = mean(diff(dist));
% clac shear
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
plot(bathy_GEBCO.dist_xsec,-bathy_GEBCO.z_xsec,'k')
axis ij tight
colorbar
colormap(redblue(25))
caxis([-1 1]*max(abs(davdx(:)))*0.5)
xlabel('Distance Along Cross-Section [km]')
ylabel('Depth [m]')
title('Vitiaz Along Channel Horizontal Shear [m/s km]')

subplot(212)
pcolor(dist_shear,sadcp.z,dxvdx)
shading interp
hold on
scatter(x_mgrid,zeros(size(x_mgrid)),100,'vk','fill')
scatter(x_mgrid,mbot,100,'^k','fill')
plot(bathy_GEBCO.dist_xsec,-bathy_GEBCO.z_xsec,'k')
axis ij tight
colorbar
colormap(redblue(25))
caxis([-1 1]*max(abs(dxvdx(:)))*0.25)
xlabel('Distance Along Cross-Section [km]')
ylabel('Depth [m]')
title('Vitiaz Across Channel Horizontal Shear [m/s km]')

print([figout 'Xsection_SADCPhorizontalShear.png'],'-dpng')


%% Plot locations of grid lines compared to the mooring line
figure('position',[0 0 800 800])
m_proj('lambert','long',[147.333 148.333],'lat',[-6.5 -5.5])
m_contour(bathy_GEBCO.lon,bathy_GEBCO.lat,bathy_GEBCO.z,-1500:100:0)
hold on
m_contour(bathy_GEBCO.lon,bathy_GEBCO.lat,bathy_GEBCO.z,[0 0],'k',...
  'LineWidth',2)
m_plot(x,y,'r')
m_scatter(sadcp.lon,sadcp.lat,10,'k','fill')
m_scatter(mlon,mlat,300,'rp','fill')
m_grid('box','fancy','tickdir','in');
title('Vitiaz Strait')

print([figout 'Xsection_map.png'],'-dpng')


%% Save bathy cross-section
% Extract bathy at coarse, grid resolution
xq = linspace(x(1),x(2),n_gp);
yq = polyval(p,xq);

% GEBCO
bathy.z = -1*interp2(bathy_GEBCO.lon,bathy_GEBCO.lat,...
  bathy_GEBCO.z,xq,yq);
bathy.SID = interp2(bathy_GEBCO.lon,bathy_GEBCO.lat,...
  bathy_GEBCO.SID,xq,yq,'nearest');
bathy.lon = xq;
bathy.lat = yq;
bathy.dist = [0 cumsum(m_lldist(xq,yq)')];

% save
save([pathin 'Xsection_Bathy.mat'],'bathy')
