% plot T-S for SBE Microcat data 

clear all
close all

main_path = '/Users/jsprintall/Documents/JanetsDocuments/southernocean/AntarcticDipoleMode/Mooring/';

%%

dpath = 'Data/SBE/ADP_'

SBE37{1} = '5821';
SBE37{2} = '6004';
SBE37{3} = '6005';

Smin = [34 33.7 34.4];
Smax = [34.7 34.7 34.7];
Tmin = [-0.5 -2 1];
Tmax = [2 2 2.1];

% load data

for ii = 1: length(SBE37)
    
    fname = ([main_path dpath SBE37{ii} '_qc.mat'])
    eval(['load ' fname])
    
    % data is in structure t
    
    figure(ii)
    plot(t.salinity,t.temperature,'r.')
    axis([Smin(ii) Smax(ii) Tmin(ii) Tmax(ii)])
    xlabel('Salinity')
    ylabel('Temperature')
    title([t.serialNo ' Mean z= ' num2str(nanmean(t.pressure))])
    
    ppath = 'Figures/';
    eval(['print -djpeg ' main_path ppath t.serialNo '_TS.jpg'])
    
end
