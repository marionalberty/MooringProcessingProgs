function moorPlotCoverage(TT,UV,time,params)
%MOORPLOTCOVERAGE Plots temporal and spatial coverage of mooring data.
% moorPlotCoverage(TT,UV,params) plots and save figures of data coverage
% for the measured temperature and velocity on the particular mooring.
% The following set of plots are generated:
% 1. A profile of the mean and standard deviation of each instrument's
% measured variable and the associated mean and standard deviation of
% observed pressure.
% 2. The full time-depth record for temperature, meridional velocity and,
% zonal velocity.
% 3.
%
% TT, SS, and UV are the data structures for temperature, salinity, and
% velocity, respectively. At this point pressure has been interpolated for
% all temperature data that did not originally have a coresponding pressure
% record.
% params contains the mooring name, file paths, and start/stop times.


%% Calculate variables needed for mean coverage profiles
% Velocity stats
u_bar = nanmean(UV.U,2);
v_bar = nanmean(UV.V,2);
p_barV = nanmean(UV.pres,2);
u_std = nanstd(UV.U,0,2);
v_std = nanstd(UV.V,0,2);
p_stdV = nanstd(UV.pres,0,2);

% Stats for velocity magnitude and direction to check for compass or
% magnitude errors
vel = UV.U + 1i*UV.V;
mag = abs(vel);
ang = angle(vel);
mag_bar = nanmean(mag,2);
ang_bar = nanmean(ang,2);
mag_std = nanstd(mag,0,2);
ang_std = nanstd(ang,0,2);

% Temperature
T_bar = nanmean(TT.temp,2);
p_barT = nanmean(TT.pres,2);
T_std = nanstd(TT.temp,0,2);
p_stdT = nanstd(TT.pres,0,2);


%% Make and print profile of mean coverage of UVT

% Find UP/DOWN ADCPs and single pt current meters
[~,~,idx] = unique(UV.snum,'stable');
unique_idx = accumarray(idx(:),(1:length(idx))',[],@(x) {sort(x)});

figure('position',[0 0 1000 600])
% Zonal velocity
ax(1) = subplot(131);
herrorbar(u_bar,p_barV,u_std,'k.')
hold on
UP = 1;
for i = 1:numel(unique_idx)
  if numel(unique_idx{i}) > 1 && UP == 1
    % UP ADCP
    errorbar(u_bar(unique_idx{i}),p_barV(unique_idx{i}),...
      p_stdV(unique_idx{i}),'rd','markerfacecolor','r')
    UP = 0;
  elseif numel(unique_idx{i}) > 1 && UP == 0
    % DOWN ADCP
    errorbar(u_bar(unique_idx{i}),p_barV(unique_idx{i}),...
      p_stdV(unique_idx{i}),'md','markerfacecolor','m')
  else
    % Single point current meter
    errorbar(u_bar(unique_idx{i}),p_barV(unique_idx{i}),...
      p_stdV(unique_idx{i}),'bd','markerfacecolor','b')
  end
end
axis ij
ylim([0 params.gridDepth])
grid on
ylabel('Pressure [dbar]')
xlabel('Zonal Velocity [m/s]')
% Meridional velocity
ax(2) = subplot(132);
herrorbar(v_bar,p_barV,v_std,'kd')
hold on
UP = 1;
for i = 1:numel(unique_idx)
  if numel(unique_idx{i}) > 1 && UP == 1
    % UP ADCP
    errorbar(v_bar(unique_idx{i}),p_barV(unique_idx{i}),...
      p_stdV(unique_idx{i}),'rd','markerfacecolor','r')
    UP = 0;
  elseif numel(unique_idx{i}) > 1 && UP == 0
    % DOWN ADCP
    errorbar(v_bar(unique_idx{i}),p_barV(unique_idx{i}),...
      p_stdV(unique_idx{i}),'md','markerfacecolor','m')
  else
    % Single point current meter
    errorbar(v_bar(unique_idx{i}),p_barV(unique_idx{i}),...
      p_stdV(unique_idx{i}),'bd','markerfacecolor','b')
  end
end
axis ij
ylim([0 params.gridDepth])
grid on
xlabel('Meridional Velocity [m/s]')
title(strrep(params.moor_name,'_',' '))
% Temperature
ax(3) = subplot(133);
herrorbar(T_bar,p_barT,T_std,'k.')
hold on
errorbar(T_bar,p_barT,p_stdT,'rd','markerfacecolor','red')
axis ij
ylim([0 params.gridDepth])
grid on
xlabel('Temperature [^\circ C]')
linkaxes(ax,'y')

% Print figure
print([params.pathfig params.channel '/' params.moor_name ...
  '/Mean_errorbars_' params.moor_name '_0'],'-dpng')


%% Make and print profile of mean coverage of magnitude and direction

figure('position',[0 0 800 600])
% Magnitude
ax(1) = subplot(121);
herrorbar(mag_bar,p_barV,mag_std,'k.')
hold on
UP = 1;
for i = 1:numel(unique_idx)
  if numel(unique_idx{i}) > 1 && UP == 1
    % UP ADCP
    errorbar(mag_bar(unique_idx{i}),p_barV(unique_idx{i}),...
      p_stdV(unique_idx{i}),'rd','markerfacecolor','r')
    UP = 0;
  elseif numel(unique_idx{i}) > 1 && UP == 0
    % DOWN ADCP
    errorbar(mag_bar(unique_idx{i}),p_barV(unique_idx{i}),...
      p_stdV(unique_idx{i}),'md','markerfacecolor','m')
  else
    % Single point current meter
    errorbar(mag_bar(unique_idx{i}),p_barV(unique_idx{i}),...
      p_stdV(unique_idx{i}),'bd','markerfacecolor','b')
  end
end
axis ij
ylim([0 params.gridDepth])
grid on
ylabel('Pressure [dbar]')
xlabel('Speed [m/s]')
title(strrep(params.moor_name,'_',' '))
% Direction
ax(2) = subplot(122);
herrorbar(ang_bar,p_barV,ang_std,'kd')
hold on
UP = 1;
for i = 1:numel(unique_idx)
  if numel(unique_idx{i}) > 1 && UP == 1
    % UP ADCP
    errorbar(ang_bar(unique_idx{i}),p_barV(unique_idx{i}),...
      p_stdV(unique_idx{i}),'rd','markerfacecolor','r')
    UP = 0;
  elseif numel(unique_idx{i}) > 1 && UP == 0
    % DOWN ADCP
    errorbar(ang_bar(unique_idx{i}),p_barV(unique_idx{i}),...
      p_stdV(unique_idx{i}),'md','markerfacecolor','m')
  else
    % Single point current meter
    errorbar(ang_bar(unique_idx{i}),p_barV(unique_idx{i}),...
      p_stdV(unique_idx{i}),'bd','markerfacecolor','b')
  end
end
axis ij
ylim([0 params.gridDepth])
grid on
xlabel('Direction [rad]')
linkaxes(ax,'y')

% Print figure
print([params.pathfig params.channel '/' params.moor_name ...
  '/Mean_errorbars_MagDir_' params.moor_name '_0'],'-dpng')


%% Make depth-time record

if ~isvector(UV.U)
  figure('position',[0 0 700 1000])
  % Zonal velocity
  ax1 = subplot(311);
  surf(repmat(time,length(UV.snum),1),UV.pres,UV.U,'edgecolor','none')
  hold on
  UP = 1;
  for i = 1:numel(unique_idx)
    if numel(unique_idx{i}) > 1 && UP == 1
      plot3(time,UV.pres(unique_idx{i}(end),:),5*ones(size(time)),'k')
      UP = 0;
    elseif numel(unique_idx{i}) > 1 && UP == 0
      plot3(time,UV.pres(unique_idx{i}(1),:),5*ones(size(time)),'k')
    end
  end
  view(2)
  colormap(ax1,redblue(25))
  datetick('x',12)
  axis ij
  cc = max(abs(UV.U(:)))*0.8;
  caxis([-cc cc])
  ylim([0 params.gridDepth])
  xlim([params.startgrid_mat params.stopgrid_mat])
  colorbar
  title([strrep(params.moor_name,'_',' ') ' Zonal (u) Velocity [m/s]'])
  % Meridional velocity
  ax2 = subplot(312);
  surf(repmat(time,length(UV.snum),1),UV.pres,UV.V,'edgecolor','none')
  hold on
  UP = 1;
  for i = 1:numel(unique_idx)
    if numel(unique_idx{i}) > 1 && UP == 1
      plot3(time,UV.pres(unique_idx{i}(end),:),5*ones(size(time)),'k')
      UP = 0;
    elseif numel(unique_idx{i}) > 1 && UP == 0
      plot3(time,UV.pres(unique_idx{i}(1),:),5*ones(size(time)),'k')
    end
  end
  view(2)
  colormap(ax2,redblue(25))
  datetick('x',12)
  axis ij
  cc = max(abs(UV.V(:)))*0.8;
  caxis([-cc cc])
  ylim([0 params.gridDepth])
  xlim([params.startgrid_mat params.stopgrid_mat])
  colorbar
  title([strrep(params.moor_name,'_',' ') ' Meridional (v) Velocity [m/s]'])
  ylabel('Pressure [dbar]')
  % Temperature
  ax3 = subplot(313);
  surf(repmat(time,length(TT.snum),1),TT.pres,TT.temp,'edgecolor','none')
  view(2)
  colormap(ax3,parula(25))
  datetick('x',12)
  axis ij
  ylim([0 params.gridDepth])
  xlim([params.startgrid_mat params.stopgrid_mat])
  colorbar
  title([strrep(params.moor_name,'_',' ') ' Temperature [^\circ C]'])
  linkaxes([ax1,ax2,ax3],'xy')
else
  figure('position',[0 0 700 330])
  % Temperature
  surf(repmat(time,length(TT.snum),1),TT.pres,TT.temp,'edgecolor','none')
  view(2)
  colormap(parula(25))
  datetick('x',12)
  axis ij
  ylim([0 params.gridDepth])
  xlim([params.startgrid_mat params.stopgrid_mat])
  colorbar
  title([strrep(params.moor_name,'_',' ') ' Temperature [^\circ C]'])
end

% Print figure
print([params.pathfig params.channel '/' params.moor_name ...
  '/' params.moor_name '_UVT_fulltimeseries'],'-dpng')


%% Make time series of individual UVT, filtered and subsampled to daily

% Filter data
dt_filt = 1;  % Desired filtered dt [days]
UV = moorFilter(UV,time,dt_filt);
TT = moorFilter(TT,time,dt_filt);

% Subsample to daily
t_start = ceil(time(1));
t_end = floor(time(end));
time_daily = t_start:t_end;
TT = moorTimeInterp(TT,time,time_daily);
UV = moorTimeInterp(UV,time,time_daily);

% Plot Temperature
nfigs = ceil(numel(TT.pdep)/8);
i_plt = 1;
for k = 1:nfigs
  n_plts = min([8 numel(TT.pdep)-i_plt+1]);
  figure('position',[0 0 500 100*n_plts])
  for i = 1:n_plts
    subplot(n_plts,1,i)
    plot(time_daily,TT.temp(i_plt,:),'k')
    grid on
    xlim([t_start t_end])
    ylim([min(TT.temp(i_plt,:)) max(TT.temp(i_plt,:))])
    datetick('x',12,'keeplimits')
    if i == 1
      title({sprintf('%s Daily Temperature [%s C] (Figure %i of %i)',...
        strrep(params.moor_name,'_',' '),'^\circ',k,nfigs);...
        sprintf('%s %s at %i m',TT.inst{i_plt},TT.snum{i_plt},...
        TT.pdep(i_plt))})
    else
      title(sprintf('%s %s at %i m',TT.inst{i_plt},TT.snum{i_plt},...
        TT.pdep(i_plt)))
    end
    i_plt = i_plt + 1;
  end
  % Print figure
  print([params.pathfig params.channel '/' params.moor_name ...
    '/' params.moor_name '_dailyTemperature_fulltimeseries_' int2str(k) ...
    'of' int2str(nfigs)],'-dpng')
end

% Plot Velocity
% Preen data
per_nan = sum(~isnan(UV.U),2)./numel(time_daily);
UV.U(per_nan < 0.1,:) = [];
UV.V(per_nan < 0.1,:) = [];
UV.pdep(per_nan < 0.1,:) = [];
UV.pres(per_nan < 0.1,:) = [];
UV.snum(per_nan < 0.1,:) = [];
UV.inst(per_nan < 0.1,:) = [];
% Find ADCP velocities and only use every 5th bin
[~,~,idx] = unique(UV.snum,'stable');
unique_idx = accumarray(idx(:),(1:length(idx))',[],@(x) {sort(x)});
plot_idx = [];
for k = 1:numel(unique_idx)
  if numel(unique_idx{k}) > 1
    plot_idx = [plot_idx; unique_idx{k}(1:5:end)];
  else
    plot_idx = [plot_idx; unique_idx{k}];
  end
end
plot_idx = sort(plot_idx);
% Set up number of figure
nfigs = ceil(numel(plot_idx)/8);
i_plt = 1;
for k = 1:nfigs
  n_plts = min([8 numel(plot_idx)-i_plt+1]);
  figure('position',[0 0 500 100*n_plts])
  for i = 1:n_plts
    subplot(n_plts,1,i)
    quiver(time_daily,zeros(size(time_daily)),UV.U(plot_idx(i_plt),:)*100,...
      UV.V(plot_idx(i_plt),:)*100,'ShowArrowHead','off','Autoscale','off')
    hold on
    quiver(t_start+20,-40,50,0,'ShowArrowHead','off','Autoscale','off')
    text(t_start+20,-80,'0.5 m/s');
    grid on
    axis([t_start t_end -100 100])
    datetick('x',12,'keeplimits')
    set(gca,'YTickLabel','')
    if i == 1
      title({sprintf('%s Daily Velocity [m/s] (Figure %i of %i)',...
        strrep(params.moor_name,'_',' '),k,nfigs);...
        sprintf('%s %s at %i m ',UV.inst{plot_idx(i_plt)},...
        UV.snum{plot_idx(i_plt)},round(UV.pdep(plot_idx(i_plt))))})
    else
      title(sprintf('%s %s at %i m',UV.inst{plot_idx(i_plt)},...
        UV.snum{plot_idx(i_plt)},round(UV.pdep(plot_idx(i_plt)))))
    end
    i_plt = i_plt +1;
  end
  % Print figure
  print([params.pathfig params.channel '/' params.moor_name ...
    '/' params.moor_name '_dailyVelocity_fulltimeseries_' int2str(k) ...
    'of' int2str(nfigs)],'-dpng')
end


%% Make UVT time-depth plot but with daily observations

time_scatU = repmat(time_daily,length(UV.snum),1);
time_scatT = repmat(time_daily,length(TT.snum),1);

if ~isvector(UV.U)
  figure('position',[0 0 700 1000])
  % Zonal velocity
  ax1 = subplot(311);
  scatter(time_scatU(:),UV.pres(:),20,UV.U(:),'filled')
  hold on
  UP = 1;
  for i = 1:numel(unique_idx)
    if numel(unique_idx{i}) > 1 && UP == 1
      plot(time_daily,UV.pres(unique_idx{i}(end),:),'k')
      UP = 0;
    elseif numel(unique_idx{i}) > 1 && UP == 0
      plot(time_daily,UV.pres(unique_idx{i}(1),:),'k')
    end
  end
  colormap(ax1,redblue(25))
  axis ij
  cc = max(abs(UV.U(:)))*0.8;
  caxis([-cc cc])
  axis([t_start t_end 0 params.gridDepth])
  datetick('x',12,'keeplimits')
  colorbar
  title([strrep(params.moor_name,'_',' ') ...
    ' Daily Zonal (u) Velocity [m/s]'])
  % Meridional velocity
  ax2 = subplot(312);
  scatter(time_scatU(:),UV.pres(:),20,UV.V(:),'filled')
  hold on
  UP = 1;
  for i = 1:numel(unique_idx)
    if numel(unique_idx{i}) > 1 && UP == 1
      plot(time_daily,UV.pres(unique_idx{i}(end),:),'k')
      UP = 0;
    elseif numel(unique_idx{i}) > 1 && UP == 0
      plot(time_daily,UV.pres(unique_idx{i}(1),:),'k')
    end
  end
  colormap(ax2,redblue(25))
  axis ij
  cc = max(abs(UV.V(:)))*0.8;
  caxis([-cc cc])
  axis([t_start t_end 0 params.gridDepth])
  datetick('x',12,'keeplimits')
  colorbar
  title([strrep(params.moor_name,'_',' ') ...
    ' Daily Meridional (v) Velocity [m/s]'])
  ylabel('Pressure [dbar]')
  % Temperature
  ax3 = subplot(313);
  scatter(time_scatT(:),TT.pres(:),20,TT.temp(:),'filled')
  colormap(ax3,parula(25))
  axis ij
  axis([t_start t_end 0 params.gridDepth])
  datetick('x',12,'keeplimits')
  colorbar
  title([strrep(params.moor_name,'_',' ') ' Daily Temperature [^\circ C]'])
  linkaxes([ax1,ax2,ax3],'xy')
else
  figure('position',[0 0 700 330])
  % Temperature
  scatter(time_scatT(:),TT.pres(:),20,TT.temp(:),'filled')
  colormap(parula(25))
  axis ij
  axis([t_start t_end 0 params.gridDepth])
  datetick('x',12,'keeplimits')
  colorbar
  title([strrep(params.moor_name,'_',' ') ' Daily Temperature [^\circ C]'])
end

% Print figure
print([params.pathfig params.channel '/' params.moor_name ...
  '/' params.moor_name '_UVT_dailytimeseries'],'-dpng')


%% Calculate variables needed for mean coverage profiles with daily values
% Velocity stats
u_bar = nanmean(UV.U,2);
v_bar = nanmean(UV.V,2);
p_barV = nanmean(UV.pres,2);
u_std = nanstd(UV.U,0,2);
v_std = nanstd(UV.V,0,2);
p_stdV = nanstd(UV.pres,0,2);

% Stats for velocity magnitude and direction to check for compass or
% magnitude errors
vel = UV.U + 1i*UV.V;
mag = abs(vel);
ang = angle(vel);
mag_bar = nanmean(mag,2);
ang_bar = nanmean(ang,2);
mag_std = nanstd(mag,0,2);
ang_std = nanstd(ang,0,2);

% Temperature
T_bar = nanmean(TT.temp,2);
p_barT = nanmean(TT.pres,2);
T_std = nanstd(TT.temp,0,2);
p_stdT = nanstd(TT.pres,0,2);


%% Make and print profile of mean coverage of UVT

% Find UP/DOWN ADCPs and single pt current meters
[~,~,idx] = unique(UV.snum,'stable');
unique_idx = accumarray(idx(:),(1:length(idx))',[],@(x) {sort(x)});

figure('position',[0 0 1000 600])
% Zonal velocity
ax(1) = subplot(131);
herrorbar(u_bar,p_barV,u_std,'k.')
hold on
UP = 1;
for i = 1:numel(unique_idx)
  if numel(unique_idx{i}) > 1 && UP == 1
    % UP ADCP
    errorbar(u_bar(unique_idx{i}),p_barV(unique_idx{i}),...
      p_stdV(unique_idx{i}),'rd','markerfacecolor','r')
    UP = 0;
  elseif numel(unique_idx{i}) > 1 && UP == 0
    % DOWN ADCP
    errorbar(u_bar(unique_idx{i}),p_barV(unique_idx{i}),...
      p_stdV(unique_idx{i}),'md','markerfacecolor','m')
  else
    % Single point current meter
    errorbar(u_bar(unique_idx{i}),p_barV(unique_idx{i}),...
      p_stdV(unique_idx{i}),'bd','markerfacecolor','b')
  end
end
axis ij
ylim([0 params.gridDepth])
grid on
ylabel('Pressure [dbar]')
xlabel('Zonal Velocity [m/s]')
% Meridional velocity
ax(2) = subplot(132);
herrorbar(v_bar,p_barV,v_std,'kd')
hold on
UP = 1;
for i = 1:numel(unique_idx)
  if numel(unique_idx{i}) > 1 && UP == 1
    % UP ADCP
    errorbar(v_bar(unique_idx{i}),p_barV(unique_idx{i}),...
      p_stdV(unique_idx{i}),'rd','markerfacecolor','r')
    UP = 0;
  elseif numel(unique_idx{i}) > 1 && UP == 0
    % DOWN ADCP
    errorbar(v_bar(unique_idx{i}),p_barV(unique_idx{i}),...
      p_stdV(unique_idx{i}),'md','markerfacecolor','m')
  else
    % Single point current meter
    errorbar(v_bar(unique_idx{i}),p_barV(unique_idx{i}),...
      p_stdV(unique_idx{i}),'bd','markerfacecolor','b')
  end
end
axis ij
ylim([0 params.gridDepth])
grid on
xlabel('Meridional Velocity [m/s]')
title([strrep(params.moor_name,'_',' ') ' from daily values'])
% Temperature
ax(3) = subplot(133);
herrorbar(T_bar,p_barT,T_std,'k.')
hold on
errorbar(T_bar,p_barT,p_stdT,'rd','markerfacecolor','red')
axis ij
ylim([0 params.gridDepth])
grid on
xlabel('Temperature [^\circ C]')
linkaxes(ax,'y')

% Print figure
print([params.pathfig params.channel '/' params.moor_name ...
  '/Mean_errorbars_' params.moor_name '_fromDaily'],'-dpng')


%% Make and print profile of mean coverage of magnitude and direction

figure('position',[0 0 800 600])
% Magnitude
ax(1) = subplot(121);
herrorbar(mag_bar,p_barV,mag_std,'k.')
hold on
UP = 1;
for i = 1:numel(unique_idx)
  if numel(unique_idx{i}) > 1 && UP == 1
    % UP ADCP
    errorbar(mag_bar(unique_idx{i}),p_barV(unique_idx{i}),...
      p_stdV(unique_idx{i}),'rd','markerfacecolor','r')
    UP = 0;
  elseif numel(unique_idx{i}) > 1 && UP == 0
    % DOWN ADCP
    errorbar(mag_bar(unique_idx{i}),p_barV(unique_idx{i}),...
      p_stdV(unique_idx{i}),'md','markerfacecolor','m')
  else
    % Single point current meter
    errorbar(mag_bar(unique_idx{i}),p_barV(unique_idx{i}),...
      p_stdV(unique_idx{i}),'bd','markerfacecolor','b')
  end
end
axis ij
ylim([0 params.gridDepth])
grid on
ylabel('Pressure [dbar]')
xlabel('Speed [m/s]')
title([strrep(params.moor_name,'_',' ') ' from daily values'])
% Direction
ax(2) = subplot(122);
herrorbar(ang_bar,p_barV,ang_std,'kd')
hold on
UP = 1;
for i = 1:numel(unique_idx)
  if numel(unique_idx{i}) > 1 && UP == 1
    % UP ADCP
    errorbar(ang_bar(unique_idx{i}),p_barV(unique_idx{i}),...
      p_stdV(unique_idx{i}),'rd','markerfacecolor','r')
    UP = 0;
  elseif numel(unique_idx{i}) > 1 && UP == 0
    % DOWN ADCP
    errorbar(ang_bar(unique_idx{i}),p_barV(unique_idx{i}),...
      p_stdV(unique_idx{i}),'md','markerfacecolor','m')
  else
    % Single point current meter
    errorbar(ang_bar(unique_idx{i}),p_barV(unique_idx{i}),...
      p_stdV(unique_idx{i}),'bd','markerfacecolor','b')
  end
end
axis ij
ylim([0 params.gridDepth])
grid on
xlabel('Direction [rad]')
linkaxes(ax,'y')

% Print figure
print([params.pathfig params.channel '/' params.moor_name ...
  '/Mean_errorbars_MagDir_' params.moor_name '_fromDaily'],'-dpng')
end