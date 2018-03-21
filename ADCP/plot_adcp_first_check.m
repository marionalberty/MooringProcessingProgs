function plot_adcp_first_check(ii,params,ADCPtype,SerNo,DeployDepth,adcp)
%======================================================================
% Raw ADCP timeseries to identify start/end time limits and large
% variations over time (N.B: Mind pitch and roll variations !) 
% Plot:  
%       - pitch & roll (identify thresholds)
%       - Raw ADCP pressure (identify large variations over time)
%       - Raw ADCP temperature (identify large variations over time)
%
% Author: C. Germineaud (cyril.germineaud@legos.obs-mip.fr), April 2016
%======================================================================

% Plot raw ADCP u,v components
display('                                      ')
display('plotting raw ADCP data: first check...')
display('--------------------------------------')


% Get all times within the time period
depth = adcp.depth(ii);
pitch = adcp.pitch(ii);
roll = adcp.roll(ii);
temp = adcp.temperature(ii);
pres = adcp.pressure(ii)/1000; % in dbar
ndate = adcp.time(ii)';

figure
set(gcf,'PaperUnits','centimeters')
xSize = 30; ySize = 21;
xLeft = (30-xSize)/2; yTop = (21-ySize)/2;
set(gcf,'Position',[xLeft yTop xSize*40 ySize*40])
set(gca,'fontsize',16)

subplot(411)
set(gca,'fontsize',16)
% % Add Statistics if needed...
% % Plot pitch mean & error bars
% errorbar(ndate(1:100:end),pitch(1:100:end),...
%     nanstd(pitch(1:100:end))*ones(size(ndate(1:100:end))),'k');hold on
% plot(ndate(end),(nanmean(pitch)),'*k','linewidth',4);hold on
% % Add median
% plot(ndate(1),(nanmedian(pitch)),'ok','linewidth',4);hold on
% 
% % Add Range Low-High
% plot(ndate(1),nanmax(pitch),'r');
% plot(ndate(1),nanmin(pitch),'r');

%pitch 
plot(ndate,pitch,'g');hold on
axis([ndate(1) ndate(end) round(min(pitch)) round(max(pitch))])
% montick('x','m',get(gca,'position'))
datetick('x',12)
ylabel('Pitch (deg)')
grid on; box on  

% add title
if strcmp(params.moor_name,'Solomon_M3') && strcmp(ADCPtype,'RDI75kHz')
    title(sprintf('%s_%s_%s_%dm',...
        params.moor_name,ADCPtype,num2str(SerNo),DeployDepth),...
        'interpreter','None','fontsize',16)
else
    title(sprintf('%s_%s_%d_%dm',...
        params.moor_name,ADCPtype,SerNo,DeployDepth),...
        'interpreter','None','fontsize',16)
end

subplot(412)
set(gca,'fontsize',16)
%roll
plot(ndate,roll,'g')
axis([ndate(1) ndate(end) round(min(roll)) round(max(roll))])
% montick('x','m',get(gca,'position'))
datetick('x',12)
ylabel('Roll (deg)')
grid on; box on  
          
subplot(413)
set(gca,'fontsize',16)
if (sum(pres)>0)
    % pressure
    plot(ndate,pres,'b')
    axis([ndate(1) ndate(end) round(min(pres)) round(max(pres))])
    % montick('x','m',get(gca,'position'))
    datetick('x',12)
    ylabel('Pressure (dbar)')
    grid on; box on  
    axis ij
    
else
    if(sum(depth(:,1))>0)
        % depth
        plot(ndate,depth(:,1),'b')
        axis([ndate(1) ndate(end) round(min(depth(:,1))) round(max(depth(:,1)))])
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
set(gca,'fontsize',16)
% temperature
plot(ndate,temp,'r')
axis([ndate(1) ndate(end) round(min(temp)) round(max(temp))])
% montick('x','m',get(gca,'position'))
datetick('x',12)
ylabel('Temperature (°C)')
grid on; box on          
orient tall

% save figure
if params.print==1
    display('saving figure...')
    
    if isfield(params,'figs_raw')
        if strcmp(params.moor_name,'Solomon_M3') && strcmp(ADCPtype,'RDI75kHz')
            nfop=sprintf('%s_%s_%s_%4.4dm_raw_first_check',params.moor_name,ADCPtype,...
                num2str(SerNo),DeployDepth);
        else
            nfop=sprintf('%s_%s_%d_%4.4dm_raw_first_check',params.moor_name,ADCPtype,...
                SerNo,DeployDepth);
        end

        nfoppng=sprintf('%s/%s/%s.png',params.figs_raw,params.moor_name,nfop);
        set(gcf,'units','centimeters')
        set(gcf,'papersize',[32 20])
        set(gcf,'paperposition',[0,0,32,20])
        print(gcf,'-dpng','-r200',nfoppng);

        if isfield(params,'dirpdf')
            nfoppdf=sprintf('%s/%s/%s.pdf',params.figs_raw,params.moor_name,...
                            params.dirpdf,nfop);
            set(gcf,'units','centimeters')
            set(gcf,'papersize',[32 20])
            set(gcf,'paperposition',[0,0,32,20])
            print(gcf,'-dpdf','-r200',nfoppdf);
        end
    end
end
display('done!')

end