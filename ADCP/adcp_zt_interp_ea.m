% Program to read the raw adcp data
% and correct for mooring motion and interpolate in pressure (depth) and
% time

% this program does this for the echo amplitude - gets vertically
% interpolated velocity from the existing file to add in


close all
clear all
% path for some routines
addpath('/Users/janetsprintall/JanetsDocuments/indo/seminars/bogor_wshop_1007/matlab_tutorial/instant_nov07/matlab/',path);
addpath('/Applications/MATLAB74/toolbox/plots/',path);

% common paths

main_path = '/Users/janetsprintall/JanetsDocuments/philippines/data/moorings/adcp_data/';


%%
% 1. INITIALIZE INPUT FOR EACH MOORING
%
% These parameters probably need to be changed for each mooring and
% deployment!

% 1. South Mindoro: 
%------------------
% Set up a regular depth grid
% bottom depth
bottom = 580;
% depth interval
zint = 10;
z_grid = [0:zint:bottom];

% A: deployment1a
%----------------
% deployment time
start_time=datenum(2007,06,09,01,00,00);
end_time = datenum(2007,06,27,21,31,00);    % have last good data at 27-Jun-2007 21:30:35
moor_name = 'SouthMindoro_deploy1a';
moor_path = 'deploy1/smin/smin_1a/process/'
% adcp data file
fname = [main_path moor_path moor_name '_8998.mat'];

% % % B: deployment1b
% % % --------------
% % % for file names
% % moor_name = 'SouthMindoro_deploy1b';
% % moor_path = 'deploy1/smin/smin_1b/process/'
% % % adcp data file
% % fname = [main_path moor_path moor_name '_8998.mat'];
% % % deployment times
% % start_time=datenum(2007,06,28,01,00,00);
% % end_time = datenum(2007,12,05,05,30,00);

% get rid of _ as this is a subscript for printing titles!
mn = strrep(moor_name,'_',' ')

%% read in ADCP data file and vertically interpolate with depth
% We will do a piecewise spline interpolation 
% [u,v] assume a slab-type surface layer
% [u_nosurf,v_nosurf] assume no fill in surface layer

load(fname)

% water depth is 580 m 
% no data at 550

[ea1,ea2,ea3,ea4] = deal(nan(size(s.time,1),length(z_grid)));

for it = 1:length(s.time)
    if mod(it - 1, 500) == 0
        disp(['Interpolating scan no. in depth ', num2str(it)]);
    end
    
    % check for good data in u (and by default, this will also reveal the good
    % data in v) and pressure
% % 
% %     jj=~isnan(s.u(it,:))& ~ isnan(s.depth(it,:));
% %     xu = s.u(it,jj);
% %     xv = s.v(it,jj);
% %     xp = s.depth(it,jj);

% for ea you do not want to get rid of the top most bins as this is the
% interesting part to the biologists where the beams reflect!!

    % now the four ea values

    xea1 = s.QCvariables(1).ea(:,it)';
    xea2 = s.QCvariables(2).ea(:,it)';
    xea3 = s.QCvariables(3).ea(:,it)';
    xea4 = s.QCvariables(4).ea(:,it)';
    xpea = s.depth(it,:);
% %     % now, add on the assumption that velocity is zero at the bottom
% %     xu = [xu 0];
% %     xv = [xv 0];
% %     xp = [xp max(z_grid)];
% %     
% %     % now we need to sort the data from shallowest to deepest, just incase
% %     % the depth array is inverted!
% %     [junk,iz] = sort(xp);
% %     xp = xp(iz);
% %     xu = xu(iz);
% %     xv = xv(iz);

% %     [junk,iz] = sort(xpea);
% %     
% %     xea1 = xea1(iz);
% %     xea2 = xea2(iz);
% %     xea3 = xea3(iz);
% %     xea4 = xea4(iz);
% %     xpea = xpea(iz);

    ea1i = interp1(xpea,xea1,z_grid,'pchip',nan);
    ea2i = interp1(xpea,xea2,z_grid,'pchip',nan);
    ea3i = interp1(xpea,xea3,z_grid,'pchip',nan);
    ea4i = interp1(xpea,xea4,z_grid,'pchip',nan);
    
    ea1(it,:) = ea1i;
    ea2(it,:) = ea2i;
    ea3(it,:) = ea3i;
    ea4(it,:) = ea4i;
    
    % this section does for velocity but only want ea
% %     % now interpolate to new depth array using a slab surface layer
% %     %
% % 
% %     ui = interp1(xp,xu,z_grid,'pchip',xu(1));   %
% %     vi = interp1(xp,xv,z_grid,'pchip',xv(1));   %
% %     % save these data that have the interpolated
% %     % surface values
% %     
% %     u(it,:) = ui;
% %     v(it,:) = vi;
% % 
% %         % now remove the data form the surface to save as no_surf
% %         isub = find(ui==ui(1));
% %         ui(isub(1:end-1))=nan;
% %         vi(isub(1:end-1)) = nan;
% % 
% %         u_nosurf(it,:) = ui;
% %         v_nosurf(it,:) = vi;

end


%% plotting

% have switch for dateform in datetick for SMIN dep1a which was ~ 3 weeks
% go to dateform = 19 if series < 1 month
% assumes 0.5 hr time series, so 48*30 = 1440
dform = 12;
if(length(s.time)<=1440)
    dform = 19;
end

% Let's plot the ea at the first 2 beams

figure(4)
clf
subplot(2,1,1)
set(gcf,'defaultaxesFontWeight','bold')
set(gcf,'defaultAxesFontSize',15)
set(gcf,'DefaultAxesLineWidth',2)
pcolor(s.time,z_grid,ea1')
shading flat
colorbar
axis ij
title(['Depth-interpolated Echo Amplitude (Beam 1): ' mn],'fontsize',12)
set(gca,'fontsize',12)
set(gca,'tickdir','out')
ylabel('depth (m)')
axis([s.time(1),s.time(end),0,max(z_grid)]);
datetick('x',dform,'keeplimits')

subplot(2,1,2)
set(gcf,'defaultaxesFontWeight','bold')
set(gcf,'defaultAxesFontSize',15)
set(gcf,'DefaultAxesLineWidth',2)
pcolor(s.time,z_grid,ea2')
shading flat
colorbar
axis ij

title(['Depth-interpolated Echo Amplitude (Beam 2): ' mn],'fontsize',12)
set(gca,'fontsize',12)
set(gca,'tickdir','out')
ylabel('depth (m)')
va = axis;
axis([s.time(1),s.time(end),0,max(z_grid)]);
datetick('x',dform,'keeplimits')

eval(['print -djpeg ',[main_path moor_path moor_name], '_zinterp_ea12.jpg'])


%% Step 3: Save the final data set

% just remove the fields that are not useful!

s = rmfield(s,'u');  % change to depth interpolated
s = rmfield(s,'v');  % changed to depth interpolated
s = rmfield(s,'qcthresh');
s = rmfield(s,'temperature');

% now make ea1 -> ea4 into structure
ea_withz(:,:,1) = ea1;
ea_withz(:,:,2) = ea2;
ea_withz(:,:,3) = ea3;
ea_withz(:,:,4) = ea4;

s = setfield(s,'ea_withz',ea_withz);
s = setfield(s,'z_grid',z_grid);

readme = 'ea_withz are depth interpolated QCvariables(1:4).ea to z_grid';
s = setfield(s,'readme',readme);

adcp_file = [main_path moor_path moor_name '_ea.mat'];
eval(['save ' adcp_file ' s']);
