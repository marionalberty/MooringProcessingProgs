function channelPlot_G(ASV,z,time,bathy,params,inst_xz,RHO,TEM,SAL)
% CHANNELPLOT Plots mean and varying velocity and transport for the
% particular interpolation and extrapolation scheme used to generate the
% channel timeseries.

%% Calc variables for plotting

% mean velocity cross-sections
ASV_bar = nanmean(ASV,3);
ASV_std = nanstd(ASV,0,3);
c_a = max(abs(ASV_bar(:)));
c_as = max(abs(ASV_std(:)));

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


%% pcolor mean x-section of ASV

figure('position',[0 0 1200 400])
% Mean Along Strait Velocity
ax1 = subplot(121);
contourf(bathy.dist,z,ASV_bar,25,'linestyle','none')
hold on
plot(bathy.dist,bathy.z,'k')
scatter(inst_xz(:,1),inst_xz(:,2),20,'k','fill')
contour(bathy.dist,z,RHO_bar,1024:1032,'k','ShowText','on')
axis tight ij
colorbar
colormap(ax1,redblue(26))
caxis([-c_a c_a])
ylabel('Depth [m]')
xlabel('Distance along section [km]')
title({params.channel;'Mean Along Strait Velocity [m/s]'})

% STD of Along Strait Velocity
ax2 = subplot(122);
contourf(bathy.dist,z,ASV_std,25,'linestyle','none')
hold on
plot(bathy.dist,bathy.z,'k')
scatter(inst_xz(:,1),inst_xz(:,2),20,'k','fill')
contour(bathy.dist,z,SAL_bar,'k','ShowText','on')
axis tight ij
colorbar
colormap(ax2,paruly(25))
caxis([0 c_as])
title('Std of Along Strait Velocity [m/s]')

print([params.pathfig params.channel '/channelGrid/meanVelXsection_' ...
  params.channel '_' params.xpassage_method 'Xpassage_'...
  params.interp_method 'Interp_' params.bottom_method '_' ...
  params.surface_method],'-dpng')


%% Mean profile of transport and std

figure('position',[0 0 600 600])
herrorbar(Tbar_z,z,Tstd_z,'b')
hold on
plot(Tbar_z,z,'b')
axis ij tight
grid on
ylabel('Depth [m]')
xlabel('Transport per unit depth [Sv/m]')
title([params.channel ', Mean Deployment Transport ' Tbar ' Sv'])

print([params.pathfig params.channel '/channelGrid/meanTransportProfile_'...
  params.channel '_' params.xpassage_method 'Xpassage_'...
  params.interp_method 'Interp_' params.bottom_method '_' ...
  params.surface_method],'-dpng')


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
  contour(time,z,RHO_perZ,1024:1032,'k','ShowText','on')
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
  params.xpassage_method 'Xpassage_' params.interp_method 'Interp_' ...
  params.bottom_method '_' params.surface_method],'-dpng')