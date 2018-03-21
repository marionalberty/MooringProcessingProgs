function plot_adcp_level0(params,s)
%======================================================================
% ADCP level 0 data timeseries:
% Plot:
%       - pitch & roll (identify thresholds)
%       - pressure/depth (identify large variations over time)
%       - temperature (identify large variations over time)
%
% Author: C. Germineaud (cyril.germineaud@legos.obs-mip.fr), July 2016
%======================================================================

set(0,'defaultaxesfontsize',14,'defaultaxeslinewidth',0.7,...
  'defaultlinelinewidth',1,'defaultpatchlinewidth',0.7,...
  'defaultFigureColor','white')

% Plot level 0 ADCP u,v components
display('                            ')
display('plotting level 0 ADCP data: ')
display('----------------------------')

figure
subplot(411)

%pitch
plot(s.time,s.pitch,'g');hold on
axis([s.time(1) s.time(end) round(min(s.pitch)) round(max(s.pitch))])
% montick('x','m',get(gca,'position'))
datetick('x',12)
ylabel('Pitch (deg)')
grid on; box on


% add title
if strcmp(s.mooringName,'Solomon_M3') && strcmp(s.meterType,'RDI75kHz')
  title(sprintf('%s_%s_%s_%dm',...
    s.mooringName,s.meterType,s.serialNo,s.plannedMeterDepth),...
    'interpreter','None')
else
  title(sprintf('%s_%s_%d_%dm',...
    s.mooringName,s.meterType,s.serialNo,s.plannedMeterDepth),...
    'interpreter','None')
end

subplot(412)
%roll
plot(s.time,s.pitch,'g')
axis([s.time(1) s.time(end) round(min(s.pitch)) round(max(s.pitch))])
% montick('x','m',get(gca,'position'))
datetick('x',12)
ylabel('Roll (deg)')
grid on; box on

subplot(413)
if (sum(s.pressure)>0)
  % pressure
  plot(s.time,s.pressure,'b')
  axis([s.time(1) s.time(end) round(min(s.pressure)) round(max(s.pressure))])
  % montick('x','m',get(gca,'position'))
  datetick('x',12)
  ylabel('Pressure (dbar)')
  grid on; box on
  axis ij
  
else
  if(sum(s.depth(1,:))>0)
    % depth
    plot(s.time,s.depth(1,:),'b')
    axis([s.time(1) s.time(end) round(min(s.depth(1,:))) ...
      round(max(s.depth(1,:)))])
    % montick('x','m',get(gca,'position'))
    datetick('x',12)
    ylabel('Depth (m)')
    grid on; box on
    axis ij
  else
    disp('no pressure/depth plot !')
  end
  
end

subplot(414)
% temperature
plot(s.time,s.temperature,'r')
axis([s.time(1) s.time(end) round(min(s.temperature)) ...
  round(max(s.temperature))])
% montick('x','m',get(gca,'position'))
datetick('x',12)
ylabel('Temperature (°C)')
grid on; box on
orient tall

% save figure
if params.print==1
  display('saving figure...')
  
  if isfield(params,'figs_lev0')
    if strcmp(s.mooringName,'Solomon_M3') && strcmp(s.meterType,'RDI75kHz')
      nfop=sprintf('%s_%s_%s_%4.4dm_lev0',s.mooringName,s.meterType,...
        s.serialNo,s.plannedMeterDepth);
    else
      nfop=sprintf('%s_%s_%d_%4.4dm_lev0',s.mooringName,s.meterType,...
        s.serialNo,s.plannedMeterDepth);
    end
    
    nfoppng=sprintf('%s/%s/%s.png',params.figs_raw,s.mooringName,nfop);
    set(gcf,'units','centimeters')
    set(gcf,'papersize',[32 20])
    set(gcf,'paperposition',[0,0,32,20])
    print(gcf,'-dpng','-r200',nfoppng);
    
    if isfield(params,'dirpdf')
      nfoppdf=sprintf('%s/%s/%s.pdf',params.figs_raw,s.mooringName,...
        nfop);
      set(gcf,'units','centimeters')
      set(gcf,'papersize',[32 20])
      set(gcf,'paperposition',[0,0,32,20])
      print(gcf,'-dpdf','-r200',nfoppdf);
    end
  end
end
display('done!')

end