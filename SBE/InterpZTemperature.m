% program to interpolate SBE data to regular grid



clear all
close all

%addpath('/Applications/MATLAB74/toolbox/plots')

main_path = '/Users/jsprintall/Documents/JanetsDocuments/southernocean/AntarcticDipoleMode/Mooring/';

        moor_name = 'ADP_';

        % path for prints and other I/O in process subdirectory
        moor_path = 'Data/SBE/';


T_file = [main_path moor_path moor_name '_allSBE.mat'];
load([T_file]);

zgrid = 140:10:250;

Tgrid = nan(length(t.time),length(zgrid));

% now set up linear interpolation and don;t forget nans'

for it = 1:length(t.time)
    if mod(it - 1, 500) == 0
        disp(['Interpolating scan no. in depth ', num2str(it)]);
    end
    
    tt = t.temperature(it,:);
    tp = t.pressure(it,:);
    if(~isnan(tt))
        ti = interp1(tp,tt,zgrid,'linear','extrap');
        Tgrid(it,:) = ti;
    end
end

time = t.time;

display('saving file')
T_file = [main_path moor_path moor_name 'SBETinterpZ.mat'];
eval(['save ' T_file ' time zgrid Tgrid']);


%%

        figure(1)
        set(gcf,'defaultaxesFontWeight','bold')
set(gcf,'defaultAxesFontSize',12)
set(gcf,'DefaultAxesLineWidth',2)
set(gcf,'DefaultLineLineWidth',2)

contourf(t.time,zgrid,Tgrid')
axis ij
colorbar
hold on
montick('x','m',get(gca,'position'))
contour(t.time,zgrid,Tgrid',[-2:1:2])
ylabel('Depth (m)')

title(['ADP Mooring Interpolated SBE Temperature with Depth'])

    outfile = [main_path 'Figures/' moor_name 'Tzinterp.jpg']
    eval(['print -djpeg ' outfile]);


        