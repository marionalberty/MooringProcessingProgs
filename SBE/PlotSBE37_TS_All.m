% plot T-S for SBE Microcat data 

clear all
close all

main_path = '/Users/jsprintall/Documents/JanetsDocuments/southernocean/AntarcticDipoleMode/Mooring/';

%%

dpath = 'Data/SBE/ADP_'

SBE37{2} = '5821';
SBE37{1} = '6004';
SBE37{3} = '6005';

Smin = [34 33.7 34.4];
Smax = [34.7 34.7 34.7];
Tmin = [-0.5 -2 1];
Tmax = [2 2 2.1];

% load data
 col = ['r','g','b'];
 sha = ['o','s','*'];
 
 mina = [33.8 34 34.4];
 
     %% set up sigma grid
    S(1) = 33.5;
    S(2) = 35.5;
    T(1) = -2;
    T(2) = 2.5;

    % grid points ofr contouring
    Sg=S(1)+[0:30]/30*(S(2)-S(1));
    Tg=T(1)+[0:30]'/30*(T(2)-T(1));


    SG=sw_dens(ones(size(Tg))*Sg,Tg*ones(size(Sg)),190)-1000;
 
 figure(1)
set(gcf,'defaultaxesFontWeight','bold')
set(gcf,'defaultAxesFontSize',12)
set(gcf,'DefaultAxesLineWidth',2)
set(gcf,'DefaultLineLineWidth',2)

 figure(2)
set(gcf,'defaultaxesFontWeight','bold')
set(gcf,'defaultAxesFontSize',12)
set(gcf,'DefaultAxesLineWidth',2)
set(gcf,'DefaultLineLineWidth',2)



for ii = 1:length(SBE37)
    
    fname = ([main_path dpath SBE37{ii} '_qc.mat'])
    eval(['load ' fname]);
    isub = find(t.salinity>=mina(ii));
    % quick triage for totally unrealistic values
    
    % data is in structure t
    figure(1)
    plot(t.salinity(isub),t.temperature(isub),[col(ii) sha(ii)]);
    xlabel('Salinity')
    ylabel('Temperature')
    text(33.56,2.2-ii*0.3,[t.serialNo ' Mean z= ' num2str(nanmean(t.pressure),4)],'color',col(ii))
    axis([33.5 35 -2 2.2])
    hold on
            if(ii==1)
            [c,h]=contour(Sg,Tg,SG,[27.5:0.1:29.5],'color',[0.5 0.5 0.5]);
            clabel(c,h,'Rotation',10,'fontsize',10,'labelspacing',200,'color',[0.5 0.5 0.5]);
            clear c h
            end
            
    figure(2)
    % use a hamming window of 5 days

% sampling time is 0.5 hour

filt_len = 5*48;

dvec = datestr(t.time);
dsub = find(dvec(:,13)=='1'&dvec(:,14)=='2'); % ie. make a daily time series centered in middle of day
ndate = t.time(dsub);

        suf = nan*ones(size(t.salinity),1);
        sup = filt_ends(hamming(filt_len),t.salinity(isub));
        suf(isub) = sup;
        % now subsample data to a daily time array
        % fill in missing data

        Sii = interp1(t.time(isub),suf(isub),ndate,'linear',NaN);
        Si = Sii;
        %         Ti(:,iz) = Tii(dsub);

daily_S = Si;
daily_time = ndate;

plot(daily_time,daily_S,'color',col(ii))
axis([daily_time(1) daily_time(end) 33.9 34.7])
hold on
    text(datenum('1-jun-2011'),34.1-ii*0.05,[t.serialNo ' Mean z= ' num2str(nanmean(t.pressure),4)],'color',col(ii))


end

figure(1)
ppath = 'Figures/';
    eval(['print -djpeg ' main_path ppath 'TS_All.jpg'])
    
    figure(2)
    
    montick('x','m',get(gca,'position'))
    eval(['print -djpeg ' main_path ppath 'SalinityTimeSeries.jpg'])
