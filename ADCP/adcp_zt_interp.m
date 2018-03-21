% Program to read the processed ADCP data
% and correct for mooring motion and interpolate in pressure (depth) and
% time

% note that when you flip the matrix upside down and then use pcolor to
% contour, you lose the first bin in the plot, but not in the actual
% retained data set.

close all
clear all
% path for some routines
%addpath('/Users/janetsprintall/JanetsDocuments/indo/seminars/bogor_wshop_1007/matlab_tutorial/instant_nov07/matlab/',path);
%addpath('/Applications/MATLAB74/toolbox/plots/',path);

% common paths
main_path = '/Users/jsprintall/Documents/JanetsDocuments/swpacific/MoorSPICE/Moorings/';


% path for figures
file_fig = 'Figures/';
file_dat = 'Data/';

moor_depl = 1;

%%
% 1. INITIALIZE INPUT FOR EACH MOORING
%
% These parameters probably need to be changed for each mooring and
% deployment!
switch moor_depl
    case 1

                % 1. St Georges East 300 khZ
        % --------------
        % for file names
        moor_name = 'StGeorgesEast';
        ADCPdir = 'ADCP/';
        % WH300 khz upward at 152 m

        WaterDepth = 1433;
        DeployDepth = 152;
        SerNo = 16832;
        ADCPtype = 'ADCP_300kHz';

        % adcp orientation
        aorient = 'up';
        lat = -4-6.174/60;
        lon = 152+31.116/60;


        % deployment times
        start_time=datenum(2012,07,20,06,10,00);
        end_time = datenum(2014,03,8,04,00,00);

        % depth interval
        zint = 10;
        z_grid = [0:zint:WaterDepth];

end

% get rid of _ as this is a subscript for printing titles!
mn = strrep(moor_name,'_',' ');

%% read in stacked ADCP data file and vertically interpolate with depth
% We will do a piecewise spline interpolation
% [u,v] assume a slab-type surface layer
% [u_nosurf,v_nosurf] assume no fill in surface layer


adcp_file = [main_path moor_path moor_name '_allvel.mat'];
eval(['load ' adcp_file]);


% implement magnetic declination
[th,r] = cart2pol(s.u,s.v);
th = th+s.magneticDeclination;
[uu,vv] = pol2cart(th,r);

[u,v] = deal(nan(length(z_grid),size(s.time,2)));
[u_nosurf,v_nosurf] = deal(nan(length(z_grid),size(s.time,2)));


itime = find(s.time==datenum(2009,2,3,12,0,0));
itime = find(s.time==datenum(2008,2,25,06,30,0));

for it = itime:itime%1:length(s.time)
    if mod(it - 1, 500) == 0
        disp(['Interpolating scan no. in depth ', num2str(it)]);
    end
    % check for good data in u (and by default, this will also reveal the good
    % data in v) and pressure

    jj=~isnan(uu(:,it)) & ~ isnan(s.depth(:,it));
    xu = uu(jj,it);
    xv = vv(jj,it);
    xp = s.depth(jj,it);
    % determine if there are data at this depth
    if ~isempty(xu)
        % now, add on the assumption that velocity is zero at the WaterDepth
        %     xu = [xu' 0];
        %     xv = [xv' 0];
        %     xp = [xp' max(z_grid)];

        xu = [xu' 0];
        xv = [xv' 0];
        xp = [xp' max(z_grid)];


        % now we need to sort the data from shallowest to deepest, just incase
        % the depth array is inverted!
        [junk,iz] = sort(xp);
        xp = xp(iz);
        xu = xu(iz);
        xv = xv(iz);
        % now interpolate to new depth array using a slab surface layer
        %
        % check for data
        % note that pchip also puts that value for the missing values in
        % the WaterDepth layer too.
        
        ui = interp1(xp,xu,z_grid,'pchip',xu(1));   %
        vi = interp1(xp,xv,z_grid,'pchip',xv(1));   %
        % save these data that have the interpolated
        % surface values

        u(:,it) = ui';
        v(:,it) = vi';

        % now remove the data from the surface and below the deepest depth to save as no_surf
        isub = find(ui==ui(1));
        ui(isub(1:end-1))=nan;
        vi(isub(1:end-1)) = nan;

        % now at WaterDepth (depth above, as have added z_grid(end) to xp
        % above)
        isub = find(z_grid>xp(end-1));
        ui(isub) = nan;
        vi(isub) = nan;
        
        u_nosurf(:,it) = ui';
        v_nosurf(:,it) = vi';
        
figure(1)
set(gcf,'defaultaxesFontWeight','bold')
set(gcf,'defaultAxesFontSize',12)
set(gcf,'DefaultAxesLineWidth',2)
set(gcf,'DefaultLineLineWidth',2)

subplot(2,2,1)
plot(v(:,it),z_grid,'k-')
hold on
vl = interp1(xp,xv,z_grid,'linear','extrap');
plot(vl(1:6),z_grid(1:6),'g-')
hold on
plot(vi,z_grid)
hold on
axis ij
plot([0 0],[0 z_grid(end)],'k:')
axis([-0.5 0.3 0 WaterDepth])
%axis([-1 0.4 0 WaterDepth])
plot(xv(1:end-1),xp(1:end-1),'rx','markersize',8)
%plot(0.2110,0,'cx','markersize',12)
xlabel('Meridional Velocity (m/s)')
ylabel('Depth')
title(['Tablas Strait ' datestr(s.time(itime))])

eval(['print -djpeg ',[main_path moor_path moor_name], '_ProfileEg.jpg'])
%shite
        
    end  % if not empty data

end
%% % Step 2. Filtering the data
% The interpolated velocity record you have just created still contains the
% tidal and inertial frequencies, but we are more interested in the
% large-scale flow.
% As a last step, let's filter the data with a 4-day filter to keep only
% the sub-inertial flow and then sub-sample to a daily record

% set up a daily time array and find indicies
% centered at mid-day from first day

dvec = datevec(s.time);
dsub = find(dvec(:,4)==12&dvec(:,5)==0);

ndate = s.time(dsub);

% Remember that filt_ends will not work with nan data so need to make an
% array of non-nan data
% nor will it work with a time series shorter than the filter!

[nz, nt] = size(u_nosurf);
[ui,vi] = deal(nan(nz,length(dsub)));
% use a hamming window of 4 days (s.time is in hours)

filt_len = 4*24;

for iz = 1:nz
    ig = find(~isnan(u_nosurf(iz,:)));
    % put limit of at least 50% of observations needed for filt
    if(length(ig)/size(u,2) >= 0.7 && nt>filt_len )
        suf = nan*ones(size(u_nosurf(iz,:)));
        svf = nan*ones(size(v_nosurf(iz,:)));
        sup = filt_ends(hamming(filt_len),u_nosurf(iz,ig));
        svp = filt_ends(hamming(filt_len),v_nosurf(iz,ig));
        suf(ig) = sup;
        svf(ig) = svp;
        % now subsample data to a daily time array

        ui(iz,:) = suf(dsub);
        vi(iz,:) = svf(dsub);
    end
end

%% plotting
% have switch for dateform in datetick for SMIN dep1a which was ~ 3 weeks
% go to dateform = 19 if series < 1 month
% assumes 0.5 hr time series, so 48*30 = 1440

dform = 12;
if(length(s.time)<=1440)
    dform = 7;
end

% plot the raw data
figure(3)
clf
subplot(2,1,1)
set(gcf,'defaultaxesFontWeight','bold')
set(gcf,'defaultAxesFontSize',15)
set(gcf,'DefaultAxesLineWidth',2)
contourf(s.time,z_grid,u_nosurf)
axis([s.time(1),s.time(end),0,max(z_grid)]);
shading flat
caxis([-1 1])
colormap('bluewhitered')
colorbar
axis ij
title(['Zonal Velocity: ' mn],'fontsize',12)
set(gca,'fontsize',12)
set(gca,'tickdir','out')
ylabel('depth (m)')
%hold on
%[c,h]=contour(s.time,z_grid,u_nosurf',[-1.5:0.5:-0.5 0.5:0.5:1.5],'color','y','linewidth',2);
%clabel(c,h);
%va = axis;
datetick('x',dform,'keeplimits')

subplot(2,1,2)
set(gcf,'defaultaxesFontWeight','bold')
set(gcf,'defaultAxesFontSize',15)
set(gcf,'DefaultAxesLineWidth',2)
contourf(s.time,z_grid,v_nosurf)
shading flat
axis ij

caxis([-1 1])
colormap('bluewhitered')
hold on
%[c,h]=contour(s.time,z_grid,v_nosurf',[-1:0.25:-0.25 0.25:0.25:1],'color','y','linewidth',2);
%clabel(c,h);
title(['Meridional Velocity: ' mn],'fontsize',12)
set(gca,'fontsize',12)
set(gca,'tickdir','out')
ylabel('depth (m)')
%va = axis;
axis([s.time(1),s.time(end),0,max(z_grid)]);
datetick('x',dform,'keeplimits')
xlabel('June 2007')
colorbar

eval(['print -djpeg ',[main_path moor_path moor_name], '_adcp.jpg'])

% Let's plot the sub-inertial zonal velocity
figure(4)
clf
subplot(2,1,1)
set(gcf,'defaultaxesFontWeight','bold')
set(gcf,'defaultAxesFontSize',15)
set(gcf,'DefaultAxesLineWidth',2)
contourf(ndate,z_grid,ui)
shading flat
caxis([-1 1])
colormap('bluewhitered')
colorbar
axis ij
title(['Zonal Velocity (4-day Hanning Filter): ' mn],'fontsize',12)
set(gca,'fontsize',12)
set(gca,'tickdir','out')
ylabel('depth (m)')
hold on
[c,h]=contour(ndate,z_grid,ui,[-1:0.25:-0.25 0.25:0.25:1],'color','y','linewidth',2);
clabel(c,h);
%va = axis;
axis([ndate(1),ndate(end),0,max(z_grid)]);
datetick('x',dform,'keeplimits')

subplot(2,1,2)
set(gcf,'defaultaxesFontWeight','bold')
set(gcf,'defaultAxesFontSize',15)
set(gcf,'DefaultAxesLineWidth',2)
contourf(ndate,z_grid,vi)
shading flat
caxis([-1 1])
colormap('bluewhitered')
hold on
[c,h]=contour(ndate,z_grid,vi,[-1:0.25:-0.25 0.25:0.25:1],'color','y','linewidth',2);
clabel(c,h);

title(['Meridional Velocity (4-day Hanning Filter): ' mn],'fontsize',12)
set(gca,'fontsize',12)
set(gca,'tickdir','out')
ylabel('depth (m)')
va = axis;
axis([ndate(1),ndate(end),0,max(z_grid)]);
axis ij
datetick('x',dform,'keeplimits')
colorbar

eval(['print -djpeg ',[main_path moor_path moor_name], '_adcp_filt.jpg'])

%% Step 3: Save the final data set

s = setfield(s,'u_daily',ui);
s = setfield(s,'v_daily',vi);
s = setfield(s,'time_daily',ndate);
s = setfield(s,'u',u);  % change to depth interpolated
s = setfield(s,'v',v);  % changed to depth interpolated
s = setfield(s,'z_grid',z_grid);

readme = '[u,v] are depth interpolated velocity, slab surface, v=0@WaterDepth';
s = setfield(s,'readme',readme);

s = setfield(s,'u_nosurface',u_nosurf);
s = setfield(s,'v_nosurface',v_nosurf);

readme = '[u_nosurface,v_nosurface] are u,v with no surface/WaterDepth interpolation';

s = setfield(s,'readme_nosurface',readme);

adcp_file = [main_path moor_path moor_name '_zinterp_vel.mat'];
eval(['save ' adcp_file ' s']);
