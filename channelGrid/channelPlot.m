function channelPlot(ASV,XSV,z,time,bathy,params,inst_xz,RHO,TEM,SAL)
% CHANNELPLOT Plots mean and varying velocity and transport for the
% particular interpolation and extrapolation scheme used to generate the
% channel timeseries.

%% Calc variables for plotting

% mean velocity cross-sections
ASV_bar = nanmean(ASV,3);
XSV_bar = nanmean(XSV,3);
ASV_std = nanstd(ASV,0,3);
XSV_std = nanstd(XSV,0,3);
c_a = max(abs(ASV_bar(:)));
c_x = max(abs(XSV_bar(:)));
c_as = max(abs(ASV_std(:)));
c_xs = max(abs(XSV_std(:)));

RHO_bar = nanmean(RHO,3);
TEM_bar = nanmean(TEM,3);
SAL_bar = nanmean(SAL,3);

ASV_T = ASV;
ASV_T(isnan(ASV_T))= 0;

% transport
dz = mean(diff(z));
Tbar_z = nanmean(squeeze(trapz(bathy.dist*1000,ASV_T,2)),2)/(1e6);
Tstd_z = nanstd(squeeze(trapz(bathy.dist*1000,ASV_T,2)),0,2)/(1e6);
Tbar = num2str(nansum(Tbar_z)*dz,2);


%% pcolor mean x-section of ASV and XSV

figure('position',[0 0 1200 800])
% Mean Along Strait Velocity
ax1 = subplot(221);
contourf(bathy.dist,z,ASV_bar,25,'linestyle','none')
hold on
patch([bathy.dist bathy.dist(end) 0],...
  [bathy.z params.sillDepth params.sillDepth],[0.7 0.7 0.7],...
  'linestyle','none')
scatter(inst_xz(:,1),inst_xz(:,2),20,'k','fill')
contour(bathy.dist,z,RHO_bar,1024:1032,'k','ShowText','on')
axis ij
ylim([0 params.sillDepth])
xlim([0 bathy.dist(end)])
colorbar
colormap(ax1,redblue(26))
caxis([-c_a c_a])
title({params.channel;'Mean Along Strait Velocity [m/s]'})
% STD of Along Strait Velocity
ax2 = subplot(222);
contourf(bathy.dist,z,ASV_std,25,'linestyle','none')
hold on
patch([bathy.dist bathy.dist(end) 0],...
  [bathy.z params.sillDepth params.sillDepth],[0.7 0.7 0.7],...
  'linestyle','none')
scatter(inst_xz(:,1),inst_xz(:,2),20,'k','fill')
axis ij
ylim([0 params.sillDepth])
xlim([0 bathy.dist(end)])
colorbar
colormap(ax2,paruly(25))
caxis([0 c_as])
title('Std of Along Strait Velocity [m/s]')

% Mean Cross Strait Velocity
ax3 = subplot(223);
contourf(bathy.dist,z,XSV_bar,25,'linestyle','none')
hold on
patch([bathy.dist bathy.dist(end) 0],...
  [bathy.z params.sillDepth params.sillDepth],[0.7 0.7 0.7],...
  'linestyle','none')
scatter(inst_xz(:,1),inst_xz(:,2),20,'k','fill')
contour(bathy.dist,z,TEM_bar,'k','ShowText','on')
axis ij
ylim([0 params.sillDepth])
xlim([0 bathy.dist(end)])
colorbar
colormap(ax3,redblue(26))
caxis([-c_x c_x])
ylabel('Depth [m]')
xlabel('Distance along section [km]')
title('Mean Cross Strait Velocity [m/s]')
% STD of Cross Strait Velocity
ax4 = subplot(224);
contourf(bathy.dist,z,XSV_std,25,'linestyle','none')
hold on
patch([bathy.dist bathy.dist(end) 0],...
  [bathy.z params.sillDepth params.sillDepth],[0.7 0.7 0.7],...
  'linestyle','none')
scatter(inst_xz(:,1),inst_xz(:,2),20,'k','fill')
contour(bathy.dist,z,SAL_bar,'k','ShowText','on')
axis ij
ylim([0 params.sillDepth])
xlim([0 bathy.dist(end)])
colorbar
colormap(ax4,paruly(25))
caxis([0 c_xs])
title('Std of Cross Strait Velocity [m/s]')

print([params.pathfig params.channel '/channelGrid/meanVelXsection_' ...
  params.channel '_' params.xpassage_method 'Xpassage_'...
  params.bottom_method '_' params.surface_method],'-dpng')


%% Mean profile of transport and std

figure('position',[0 0 600 600])
herrorbar(Tbar_z,z,Tstd_z,'b')
hold on
plot(Tbar_z,z,'b')
axis ij
ylim([0 params.sillDepth])
grid on
ylabel('Depth [m]')
xlabel('Transport per unit depth [Sv/m]')
title([params.channel ', Mean Deployment Transport ' Tbar ' Sv'])

print([params.pathfig params.channel '/channelGrid/meanTransportProfile_'...
  params.channel '_' params.xpassage_method 'Xpassage_'...
  params.bottom_method '_' params.surface_method],'-dpng')


%% Plot Transport per unit depth

% Calc transport per unit depth
T_perZ = squeeze(trapz(bathy.dist*1000,ASV_T,2))/(1e6);
RHO_perZ = squeeze(nanmean(RHO,2));

figure('position',[0 0 1200 400])
pcolor(time,z,T_perZ)
hold on
if min(RHO_bar(:)) < 1026
  contour(time,z,RHO_perZ,[1027 1027],'k','ShowText','on')
else
  contour(time,z,RHO_perZ,1024:0.5:1032,'k','ShowText','on')
end
shading interp
clim = max(abs(T_perZ(:)));
colormap(redblue(26))
colorbar
caxis([-clim clim])
datetick('x',12,'keeplimits')
axis ij
ylim([0 params.sillDepth])
ylabel('Depth [m]')
title([params.channel ', Transport per unit depth [Sv/m]'])

print([params.pathfig params.channel '/channelGrid/'...
  'transportPerDepthTimeseries_' params.channel '_' ...
  params.xpassage_method 'Xpassage_' params.bottom_method '_' ...
  params.surface_method],'-dpng')