function moorPlotFinal(UV,SS,TT,params)
% MOORPLOTFINAL Plots the final mooring timeseries after extrapolation to
% the surface and bottom is complete. The function generates and saves the
% following figures:
% U & V in time and depth
% Mean profile of U & V with std bars
% 

%% Calculate variables need for mean profiles

% Velocity stats
u_bar = nanmean(UV.U,2);
v_bar = nanmean(UV.V,2);
u_std = nanstd(UV.U,0,2);
v_std = nanstd(UV.V,0,2);

% Density
rho = sw_dens(SS.psal,TT.temp,TT.z);


%% Mean profile of U & V

figure('position',[0 0 800 600])
% Zonal Velocity
subplot(121)
herrorbar(u_bar,UV.z,u_std,'b')
hold on
plot(u_bar,UV.z,'b')
axis ij tight
grid on
ylabel('Depth [m]')
xlabel('Zonal Velocity [m/s]')
title(strrep(params.moor_name,'_',' '))
% Meridional velocity
subplot(122)
herrorbar(v_bar,UV.z,v_std,'b')
hold on
plot(v_bar,UV.z,'b')
axis ij tight
grid on
xlabel('Meridional Velocity [m/s]')

% Print figure
print([params.pathfig params.channel '/' params.moor_name ...
  '/extrap_mooringTimeseries/MeanUV_' params.moor_name '_' ...
  params.bottom_method '_' params.surface_method],'-dpng')


%% Pcolor of U & V

figure('position',[0 0 700 660])
% Zonal Velocity
subplot(211)
pcolor(UV.time,UV.z,UV.U)
shading flat
axis ij tight
datetick('x',12,'keeplimits')
colorbar
colormap(redblue(25))
caxis([-max(abs(UV.U(:))) max(abs(UV.U(:)))])
hold on
contour(TT.time,TT.z,rho,'k')
ylabel('Depth [m]')
title({strrep(params.moor_name,'_',' ');'Zonal Velocity [m/s]'})
% Meridional velocity
subplot(212)
pcolor(UV.time,UV.z,UV.V)
shading flat
axis ij tight
datetick('x',12,'keeplimits')
colorbar
colormap(redblue(25))
caxis([-max(abs(UV.U(:))) max(abs(UV.U(:)))])
hold on
contour(TT.time,TT.z,rho,'k')
ylabel('Depth [m]')
title('Meridional Velocity [m/s]')

% Print figure
print([params.pathfig params.channel '/' params.moor_name ...
  '/extrap_mooringTimeseries/UVtimeseries_' params.moor_name '_' ...
  params.bottom_method '_' params.surface_method],'-dpng')