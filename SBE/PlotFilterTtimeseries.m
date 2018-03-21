% program tofilter then plot SBE temperature data with overlay of potential
% density



clear all
close all

%addpath('/Applications/MATLAB74/toolbox/plots')

main_path = '/Users/jsprintall/Documents/JanetsDocuments/southernocean/AntarcticDipoleMode/Mooring/';

moor_name = 'ADP_';

% path for prints and other I/O in process subdirectory
moor_path = 'Data/SBE/';


T_file = [main_path moor_path moor_name 'SBETinterpZ.mat'];
load(T_file);

%% use a hamming window of 5 days

% sampling time is 0.5 hour

filt_len = 5*48;

[nt,nz] = size(Tgrid);
dvec = datestr(time);
dsub = find(dvec(:,13)=='1'&dvec(:,14)=='2'); % ie. make a daily time series centered in middle of day

ndate = time(dsub);

for iz = 1:nz
    ig = find(~isnan(Tgrid(:,iz)));
    if(length(ig)/nt >= 0.5 & nt>filt_len )
        suf = nan*ones(size(Tgrid(:,iz)));
        sup = filt_ends(hamming(filt_len),Tgrid(ig,iz));
        suf(ig) = sup;
        % now subsample data to a daily time array
        % fill in missing data

        Tii = interp1(time(ig),suf(ig),ndate,'linear',NaN);
        Ti(:,iz) = Tii;
        %         Ti(:,iz) = Tii(dsub);

    end
end
daily_Tgrid = Ti;
daily_time = ndate;



%%

figure(1)
set(gcf,'defaultaxesFontWeight','bold')
set(gcf,'defaultAxesFontSize',12)
set(gcf,'DefaultAxesLineWidth',2)
set(gcf,'DefaultLineLineWidth',2)

contourf(daily_time,zgrid,daily_Tgrid',[-2:0.25:2])
axis ij
caxis([-2 2]);
colorbar
hold on
montick('x','m',get(gca,'position'))
[c,h]=contour(daily_time,zgrid,daily_Tgrid',[-2:1:2],'color','w','linewidth',2);
clabel(c,h,'fontsize',12,'color','w','labelspacing',300,'fontweight','bold');

ylabel('pressure')

title(['ADP Mooring Temperature'])

outfile = [main_path 'Figures/' moor_name 'Tzinterp_daily.jpg']
eval(['print -djpeg ' outfile]);


