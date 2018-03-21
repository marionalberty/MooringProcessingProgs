function plot_adcp_qc(params,nbins,s,uu,vv,qc)
%====================================================================
% Key: Plot ADCP mooring quality control tests 
%
% Author: J. Sprintall (jsprintall@ucsd.edu), Nov 2007
%
% Edited and Modified as a function by: 
% C. Germineaud (cyril.germineaud@legos.obs-mip.fr), Aug 2015
%====================================================================

% Plot raw ADCP qced data
display('                          ')
display('plotting qc data...')
display('--------------------------')

figure(103)
pcolor(s.time,1:nbins,s.u)
shading flat
%colormap('redblue')
colormap('jet')
colorbar
%caxis([-max(abs(s.u(:))) max(abs(s.u(:)))])
caxis([-2 2])
title(sprintf('Zonal Velocity QC (not mooring motion corrected) %s_%s_%s',...
              params.moor_name,s.meterType,s.serialNo),...
              'interpreter','None','fontsize',10)
          
% montick('x','m',get(gca,'position'))
datetick('x',12)
axis([s.time(1) s.time(end) 1 nbins])

if params.print==1
    if isfield(params,'figs_lev0')
        %ss=input('print figures ? ','s');
        %if ss=='y'
            nfop=sprintf('%s_%s_%s_%4.4d_%d_adcp_uqc',params.moor_name,s.meterType,...
                         s.serialNo,s.plannedMeterDepth,params.qcLevel);
            if isfield(params,'dirpng')
                nfoppng=sprintf('%s/%s.png',params.dirpng,nfop);
                set(gcf,'units','centimeters')
                set(gcf,'papersize',[32 20])
                set(gcf,'paperposition',[0,0,32,20])
                print(gcf,'-dpng','-r200',nfoppng);
            end
            if isfield(params,'dirpdf')
                nfoppdf=sprintf('%s/%s.pdf',params.dirpdf,nfop);
                set(gcf,'units','centimeters')
                set(gcf,'papersize',[32 20])
                set(gcf,'paperposition',[0,0,32,20])
                print(gcf,'-dpdf','-r200',nfoppdf);
            end
        %end
    end
end

close

figure(104)
pcolor(s.time,1:nbins,s.v)
shading flat
%colormap('redblue')
colormap('jet')
colorbar
%caxis([-max(abs(s.v(:))) max(abs(s.v(:)))])
caxis([-2 2])
% montick('x','m',get(gca,'position'))
datetick('x',12)
axis([s.time(1) s.time(end) 1 nbins])
title(sprintf('Meridional Velocity QC (not mooring motion corrected) %s_%s_%s',...
              params.moor_name,s.meterType,s.serialNo),...
              'interpreter','None','fontsize',10)

if params.print==1
    if isfield(params,'figs_lev0')
        %ss=input('print figures ? ','s');
        %if ss=='y'
        nfop=sprintf('%s_%s_%s_%4.4d_%d_adcp_vqc',params.moor_name,s.meterType,...
            s.serialNo,s.plannedMeterDepth,params.qcLevel);
        if isfield(params,'dirpng')
            nfoppng=sprintf('%s/%s.png',params.dirpng,nfop);
            set(gcf,'units','centimeters')
            set(gcf,'papersize',[32 20])
            set(gcf,'paperposition',[0,0,32,20])
            print(gcf,'-dpng','-r200',nfoppng);
        end
        if isfield(params,'dirpdf')
            nfoppdf=sprintf('%s/%s.pdf',params.dirpdf,nfop);
            set(gcf,'units','centimeters')
            set(gcf,'papersize',[32 20])
            set(gcf,'paperposition',[0,0,32,20])
            print(gcf,'-dpdf','-r200',nfoppdf);
        end
        %end
    end
end

close
    
figure(105)
pcolor(s.time,1:nbins,uu)
shading flat
%colormap('redblue')
colormap('jet')
colorbar
%caxis([-max(abs(uu(:))) max(abs(uu(:)))])
caxis([-2 2])
% montick('x','m',get(gca,'position'))
datetick('x',12)
axis([s.time(1) s.time(end) 1 nbins])
title(sprintf('Zonal Velocity Raw (not mooring motion corrected) %s_%s_%s',...
              params.moor_name,s.meterType,s.serialNo),...
              'interpreter','None','fontsize',10)

if params.print==1
    if params.print==1
        if isfield(params,'figs_lev0')
            %ss=input('print figures ? ','s');
            %if ss=='y'
            nfop=sprintf('%s_%s_%s_%4.4d_%d_adcp_uraw',params.moor_name,s.meterType,...
                s.serialNo,s.plannedMeterDepth,params.qcLevel);
            if isfield(params,'dirpng')
                nfoppng=sprintf('%s/%s.png',params.dirpng,nfop);
                set(gcf,'units','centimeters')
                set(gcf,'papersize',[32 20])
                set(gcf,'paperposition',[0,0,32,20])
                print(gcf,'-dpng','-r200',nfoppng);
            end
            if isfield(params,'dirpdf')
                nfoppdf=sprintf('%s/%s.pdf',params.dirpdf,nfop);
                set(gcf,'units','centimeters')
                set(gcf,'papersize',[32 20])
                set(gcf,'paperposition',[0,0,32,20])
                print(gcf,'-dpdf','-r200',nfoppdf);
            end
            %end
        end
    end
end

close

figure(106)
pcolor(s.time,1:nbins,vv)
shading flat
%colormap('redblue')
colormap('jet')
colorbar
%caxis([-max(abs(vv(:))) max(abs(vv(:)))])
caxis([-2 2])
% montick('x','m',get(gca,'position'))
datetick('x',12)
axis([s.time(1) s.time(end) 1 nbins])
title(sprintf('Meridional Velocity Raw (not mooring motion corrected) %s_%s_%s',...
              params.moor_name,s.meterType,s.serialNo),...
              'interpreter','None','fontsize',10)

if params.print==1
    if params.print==1
        if isfield(params,'figs_lev0')
            %ss=input('print figures ? ','s');
            %if ss=='y'
            nfop=sprintf('%s_%s_%s_%4.4d_%d_adcp_vraw',params.moor_name,s.meterType,...
                s.serialNo,s.plannedMeterDepth,params.qcLevel);
            if isfield(params,'dirpng')
                nfoppng=sprintf('%s/%s.png',params.dirpng,nfop);
                set(gcf,'units','centimeters')
                set(gcf,'papersize',[32 20])
                set(gcf,'paperposition',[0,0,32,20])
                print(gcf,'-dpng','-r200',nfoppng);
            end
            if isfield(params,'dirpdf')
                nfoppdf=sprintf('%s/%s.pdf',params.dirpdf,nfop);
                set(gcf,'units','centimeters')
                set(gcf,'papersize',[32 20])
                set(gcf,'paperposition',[0,0,32,20])
                print(gcf,'-dpdf','-r200',nfoppdf);
            end
            %end
        end
    end
end   
close

%% 4. PLOT THE ADCP QC DIAGNOSTICS

% Plot ADCP qc diagnostics
display('                          ')
display('plotting qc diagnostics...')
display('--------------------------')

figure(1)
subplot(411)
pcolor(s.time,1:nbins,(qc(1).ea));
shading flat
colorbar
% montick('x','m',get(gca,'position'))
datetick('x',12)
axis([s.time(1) s.time(end) 1 nbins])
set(gca,'tickdir','out')
ylabel('ADCP bin')
title(sprintf('Echo Amplitude Beam 1: %s_%s_%s',...
              params.moor_name,s.meterType,s.serialNo),...
              'interpreter','None','fontsize',10)
subplot(412)
pcolor(s.time,1:nbins,(qc(2).ea));
shading flat
colorbar
% montick('x','m',get(gca,'position'))
datetick('x',12)
axis([s.time(1) s.time(end) 1 nbins])
set(gca,'tickdir','out')
ylabel('ADCP bin')
title(sprintf('Echo Amplitude Beam 2: %s_%s_%s',...
              params.moor_name,s.meterType,s.serialNo),...
              'interpreter','None','fontsize',10)
subplot(413)
pcolor(s.time,1:nbins,(qc(3).ea));
shading flat
colorbar
% montick('x','m',get(gca,'position'))
datetick('x',12)
axis([s.time(1) s.time(end) 1 nbins])
set(gca,'tickdir','out')
ylabel('ADCP bin')
title(sprintf('Echo Amplitude Beam 3: %s_%s_%s',...
              params.moor_name,s.meterType,s.serialNo),...
              'interpreter','None','fontsize',10)
subplot(414)
pcolor(s.time,1:nbins,squeeze(qc(4).ea));
shading flat
colorbar
% montick('x','m',get(gca,'position'))
datetick('x',12)
axis([s.time(1) s.time(end) 1 nbins])
set(gca,'tickdir','out')
ylabel('ADCP bin')
title(sprintf('Echo Amplitude Beam 4: %s_%s_%s',...
              params.moor_name,s.meterType,s.serialNo),...
              'interpreter','None','fontsize',10)
orient tall

if params.print==1
    if params.print==1
        if isfield(params,'figs_lev0')
            %ss=input('print figures ? ','s');
            %if ss=='y'
            nfop=sprintf('%s_%s_%s_%4.4d_%d_echoamp',params.moor_name,s.meterType,...
                s.serialNo,s.plannedMeterDepth,params.qcLevel);
            if isfield(params,'dirpng')
                nfoppng=sprintf('%s/%s.png',params.dirpng,nfop);
                set(gcf,'units','centimeters')
                set(gcf,'papersize',[32 20])
                set(gcf,'paperposition',[0,0,32,20])
                print(gcf,'-dpng','-r200',nfoppng);
            end
            if isfield(params,'dirpdf')
                nfoppdf=sprintf('%s/%s.pdf',params.dirpdf,nfop);
                set(gcf,'units','centimeters')
                set(gcf,'papersize',[32 20])
                set(gcf,'paperposition',[0,0,32,20])
                print(gcf,'-dpdf','-r200',nfoppdf);
            end
            %end
        end
    end
end   

close

% %  In the other coordinate frames (ADCP, Ship and Earth Coordinates),
% the four Percent Good values represent (in order):
% 1) The percentage of good three beam solutions (one beam rejected);
% 2) The percentage of good transformations (error velocity threshold not exceeded);
% 3) The percentage of measurements where more than one beam was bad; and
% 4) The percentage of measurements with four beam solutions.

%1 adcp percent good
figure(2)
subplot(411)
pcolor(s.time,1:nbins,qc(1).pg);
shading flat
colorbar
% montick('x','m',get(gca,'position'))
datetick('x',12)
axis([s.time(1) s.time(end) 1 nbins])
set(gca,'tickdir','out')
ylabel('ADCP bin')
title(sprintf('Percent Good Beam 1: good of 3-beam solutions with 1-beam rejected: %s_%s_%s',...
              params.moor_name,s.meterType,s.serialNo),...
              'interpreter','None','fontsize',10)
subplot(412)
pcolor(s.time,1:nbins,qc(2).pg);
shading flat
colorbar
% montick('x','m',get(gca,'position'))
datetick('x',12)
axis([s.time(1) s.time(end) 1 nbins])
set(gca,'tickdir','out')
ylabel('ADCP bin')
title(sprintf('Percent Good Beam 2: Good transformation where ErVel not exceeded: %s_%s_%s',...
              params.moor_name,s.meterType,s.serialNo),...
              'interpreter','None','fontsize',10)

% 2) The percentage of good transformations (error velocity threshold not exceeded);
% 3) The percentage of measurements where more than one beam was bad; and
% 4) The percentage of measurements with four beam solutions.
subplot(413)
pcolor(s.time,1:nbins,qc(3).pg);
shading flat
colorbar
% montick('x','m',get(gca,'position'))
datetick('x',12)
axis([s.time(1) s.time(end) 1 nbins])
set(gca,'tickdir','out')
ylabel('ADCP bin')
title(sprintf('Percent Good Beam 3: where more than one beam bad: %s_%s_%s',...
              params.moor_name,s.meterType,s.serialNo),...
              'interpreter','None','fontsize',10)
subplot(414)
pcolor(s.time,1:nbins,qc(4).pg);
shading flat
colorbar
% montick('x','m',get(gca,'position'))
datetick('x',12)
axis([s.time(1) s.time(end) 1 nbins])
set(gca,'tickdir','out')
ylabel('ADCP bin')
title(sprintf('Percent Good Beam 4: with good 4-beam solutions: %s_%s_%s',...
              params.moor_name,s.meterType,s.serialNo),...
              'interpreter','None','fontsize',10)
orient tall

if params.print==1
    if params.print==1
        if isfield(params,'figs_lev0')
            %ss=input('print figures ? ','s');
            %if ss=='y'
            nfop=sprintf('%s_%s_%s_%4.4d_%d_percgood',params.moor_name,s.meterType,...
                s.serialNo,s.plannedMeterDepth,params.qcLevel);
            if isfield(params,'dirpng')
                nfoppng=sprintf('%s/%s.png',params.dirpng,nfop);
                set(gcf,'units','centimeters')
                set(gcf,'papersize',[32 20])
                set(gcf,'paperposition',[0,0,32,20])
                print(gcf,'-dpng','-r200',nfoppng);
            end
            if isfield(params,'dirpdf')
                nfoppdf=sprintf('%s/%s.pdf',params.dirpdf,nfop);
                set(gcf,'units','centimeters')
                set(gcf,'papersize',[32 20])
                set(gcf,'paperposition',[0,0,32,20])
                print(gcf,'-dpdf','-r200',nfoppdf);
            end
            %end
        end
    end
end   
close

% pitch and roll
figure(4)
subplot(411)
plot(s.time,s.pitch)
% montick('x','m',get(gca,'position'))
datetick('x',12)
ylabel('degrees')
axis([s.time(1) s.time(end) min(s.pitch) max(s.pitch)])
title(sprintf('ADCP Pitch: %s_%s_%s',...
              params.moor_name,s.meterType,s.serialNo),...
              'interpreter','None','fontsize',10)

subplot(412)
plot(s.time,s.pitch_std)
% montick('x','m',get(gca,'position'))
datetick('x',12)
axis([s.time(1) s.time(end) min(s.pitch_std) max(s.pitch_std)])
ylabel('degrees')
title(sprintf('ADCP Pitch Standard deviation: %s_%s_%s',...
              params.moor_name,s.meterType,s.serialNo),...
              'interpreter','None','fontsize',10)

subplot(413)
plot(s.time,s.roll)
axis([s.time(1) s.time(end) min(s.roll) max(s.roll)])
% montick('x','m',get(gca,'position'))
datetick('x',12)
ylabel('degrees')
title(sprintf('ADCP Roll: %s_%s_%s',...
              params.moor_name,s.meterType,s.serialNo),...
              'interpreter','None','fontsize',10)
subplot(414)
plot(s.time,s.roll_std)
axis([s.time(1) s.time(end) min(s.roll_std) max(s.roll_std)])
% montick('x','m',get(gca,'position'))
datetick('x',12)
ylabel('degrees')
title(sprintf('ADCP Roll Standard Deviation %s_%s_%s',...
              params.moor_name,s.meterType,s.serialNo),...
              'interpreter','None','fontsize',10)
orient tall

if params.print==1
    if params.print==1
        if isfield(params,'figs_lev0')
            %ss=input('print figures ? ','s');
            %if ss=='y'
            nfop=sprintf('%s_%s_%s_%4.4d_%d_pitchroll',params.moor_name,s.meterType,...
                s.serialNo,s.plannedMeterDepth,params.qcLevel);
            if isfield(params,'dirpng')
                nfoppng=sprintf('%s/%s.png',params.dirpng,nfop);
                set(gcf,'units','centimeters')
                set(gcf,'papersize',[32 20])
                set(gcf,'paperposition',[0,0,32,20])
                print(gcf,'-dpng','-r200',nfoppng);
            end
            if isfield(params,'dirpdf')
                nfoppdf=sprintf('%s/%s.pdf',params.dirpdf,nfop);
                set(gcf,'units','centimeters')
                set(gcf,'papersize',[32 20])
                set(gcf,'paperposition',[0,0,32,20])
                print(gcf,'-dpdf','-r200',nfoppdf);
            end
            %end
        end
    end
end   
close

figure(61)
subplot(211)
plot(s.time,s.heading)
axis([s.time(1) s.time(end) min(s.heading) max(s.heading)])
% montick('x','m',get(gca,'position'))
datetick('x',12)
ylabel('degrees')
title(sprintf('ADCP Heading: %s_%s_%s',...
              params.moor_name,s.meterType,s.serialNo),...
              'interpreter','None','fontsize',10)

subplot(212)
plot(s.time,s.heading_std)
axis([s.time(1) s.time(end) min(s.heading_std) max(s.heading_std)])
% montick('x','m',get(gca,'position'))
datetick('x',12)
ylabel('degrees')
title(sprintf('ADCP Heading Standard Deviation %s_%s_%s',...
              params.moor_name,s.meterType,s.serialNo),...
              'interpreter','None','fontsize',10)
orient tall

if params.print==1
    if params.print==1
        if isfield(params,'figs_lev0')
            %ss=input('print figures ? ','s');
            %if ss=='y'
            nfop=sprintf('%s_%s_%s_%4.4d_%d_heading',params.moor_name,s.meterType,...
                s.serialNo,s.plannedMeterDepth,params.qcLevel);
            if isfield(params,'dirpng')
                nfoppng=sprintf('%s/%s.png',params.dirpng,nfop);
                set(gcf,'units','centimeters')
                set(gcf,'papersize',[32 20])
                set(gcf,'paperposition',[0,0,32,20])
                print(gcf,'-dpng','-r200',nfoppng);
            end
            if isfield(params,'dirpdf')
                nfoppdf=sprintf('%s/%s.pdf',params.dirpdf,nfop);
                set(gcf,'units','centimeters')
                set(gcf,'papersize',[32 20])
                set(gcf,'paperposition',[0,0,32,20])
                print(gcf,'-dpdf','-r200',nfoppdf);
            end
            %end
        end
    end
end 
close

figure(6)
plot(s.time,s.temperature,'r')
axis([s.time(1) s.time(end) min(s.temperature) max(s.temperature)])
% montick('x','m',get(gca,'position'))
datetick('x',12)
ylabel('Temperature (°C)')
title(sprintf('ADCP Temperature: %s_%s_%s',...
              params.moor_name,s.meterType,s.serialNo),...
              'interpreter','None','fontsize',10)

if params.print==1
    if params.print==1
        if isfield(params,'figs_lev0')
            %ss=input('print figures ? ','s');
            %if ss=='y'
            nfop=sprintf('%s_%s_%s_%4.4d_%d_temperature',params.moor_name,s.meterType,...
                s.serialNo,s.plannedMeterDepth,params.qcLevel);
            if isfield(params,'dirpng')
                nfoppng=sprintf('%s/%s.png',params.dirpng,nfop);
                set(gcf,'units','centimeters')
                set(gcf,'papersize',[32 20])
                set(gcf,'paperposition',[0,0,32,20])
                print(gcf,'-dpng','-r200',nfoppng);
            end
            if isfield(params,'dirpdf')
                nfoppdf=sprintf('%s/%s.pdf',params.dirpdf,nfop);
                set(gcf,'units','centimeters')
                set(gcf,'papersize',[32 20])
                set(gcf,'paperposition',[0,0,32,20])
                print(gcf,'-dpdf','-r200',nfoppdf);
            end
            %end
        end
    end
end   
close

if (sum(s.pressure)>0)
    % pressure
    figure(5)
    if strcmp(params.moor_name,'Solomon_M1') ||... 
        strcmp(params.moor_name,'Solomon_M2b') ||...
        strcmp(params.moor_name,'Solomon_M3') ||...
        strcmp(params.moor_name,'StGeorgesEast') ||...
        strcmp(params.moor_name,'StGeorgesWest') ||...
        strcmp(params.moor_name,'VitiazMiddle')
        plot(s.time,s.depth(:,1),'b')
        axis([s.time(1) s.time(end) min(s.depth(:,1)) max(s.depth(:,1))])
        axis ij
        datetick('x',12)
        ylabel('Depth (m)')
        title(sprintf('ADCP Depth: %s_%s_%s',...
            params.moor_name,s.meterType,s.serialNo),...
            'interpreter','None','fontsize',10)
    else
        plot(s.time,s.pressure)
        axis([s.time(1) s.time(end) min(s.pressure) max(s.pressure)])
        % montick('x','m',get(gca,'position'))
        datetick('x',12)
        ylabel('Pressure (dbar)')
        axis ij
        title(sprintf('ADCP Pressure: %s_%s_%s',...
            params.moor_name,s.meterType,s.serialNo),...
            'interpreter','None','fontsize',10)  
    end

    if params.print==1
        if isfield(params,'figs_lev0')
            %ss=input('print figures ? ','s');
            %if ss=='y'
            if strcmp(params.moor_name,'Solomon_M1')
                nfop=sprintf('%s_%s_%s_%4.4d_%d_depth',params.moor_name,s.meterType,...
                    s.serialNo,s.plannedMeterDepth,params.qcLevel);
            else
                nfop=sprintf('%s_%s_%s_%4.4d_%d_pressure',params.moor_name,s.meterType,...
                    s.serialNo,s.plannedMeterDepth,params.qcLevel);
            end
            if isfield(params,'dirpng')
                nfoppng=sprintf('%s/%s.png',params.dirpng,nfop);
                set(gcf,'units','centimeters')
                set(gcf,'papersize',[32 20])
                set(gcf,'paperposition',[0,0,32,20])
                print(gcf,'-dpng','-r200',nfoppng);
            end
            if isfield(params,'dirpdf')
                nfoppdf=sprintf('%s/%s.pdf',params.dirpdf,nfop);
                set(gcf,'units','centimeters')
                set(gcf,'papersize',[32 20])
                set(gcf,'paperposition',[0,0,32,20])
                print(gcf,'-dpdf','-r200',nfoppdf);
            end
            %end
        end
    end

    
else
    if(sum(s.depth(:,1))>0)
        % depth
        figure(5)
        plot(s.time,s.depth(:,1))
        
        axis([s.time(1) s.time(end) min(s.depth(:,1)) max(s.depth(:,1))])
        axis ij
        % montick('x','m',get(gca,'position'))
        datetick('x',12)
        ylabel('Depth (m)')
        title(sprintf('ADCP Depth: %s_%s_%s',...
              params.moor_name,s.meterType,s.serialNo),...
              'interpreter','None','fontsize',10)
        
        if params.print==1
            if isfield(params,'figs_lev0')
                %ss=input('print figures ? ','s');
                %if ss=='y'
                nfop=sprintf('%s_%s_%s_%4.4d_%d_depth',params.moor_name,s.meterType,...
                    s.serialNo,s.plannedMeterDepth,params.qcLevel);
                if isfield(params,'dirpng')
                    nfoppng=sprintf('%s/%s.png',params.dirpng,nfop);
                    set(gcf,'units','centimeters')
                    set(gcf,'papersize',[32 20])
                    set(gcf,'paperposition',[0,0,32,20])
                    print(gcf,'-dpng','-r200',nfoppng);
                end
                if isfield(params,'dirpdf')
                    nfoppdf=sprintf('%s/%s.pdf',params.dirpdf,nfop);
                    set(gcf,'units','centimeters')
                    set(gcf,'papersize',[32 20])
                    set(gcf,'paperposition',[0,0,32,20])
                    print(gcf,'-dpdf','-r200',nfoppdf);
                end
                %end
            end
        end
        
    else
        disp('no pressure/depth plot !')
    end

end

close all

display('All figures plotted successfully !')

end

    
