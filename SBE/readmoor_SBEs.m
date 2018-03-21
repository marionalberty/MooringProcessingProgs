% read in SBEs from mooring data

clear all
close all

addpath('/Applications/MATLAB74/toolbox/plots')


main_path = '/Users/jsprintall/Documents/JanetsDocuments/swpacific/MoorSPICE/Moorings/';

moor_depl = 5;

% path for figures
file_fig = 'Figures/'
file_dat = 'Data/'
% path for input file from raw subdirectory

%%
% 1. INITIALIZE INPUT FOR EACH MOORING
%
% These parameters probably need to be changed for each mooring and
% deployment!

switch moor_depl
    case 1
        % ST Georges West
        moor_name = 'StGeorgesWest';
        mn = 'SGW'


        lat = -4-6.174/60;
        lon = 152+33.804/60;
        waterz = 1433;

        % deployment times
        start_time=datenum(2012,07,20,06,00,00);
        end_time = datenum(2014,03,7,21,00,00);

        % need to get rid of headers on SBE56 csv files except date,time,temp

        inst{1} = 'SBE56';
        fname{1} = 'SGW-0600m_SBE56-01024_2014-03-09.asc.csv'
        sno{1} = '01024';
        plannedz{1} = 600;


        inst{2} = 'SBE39';
        fname{2} = 'SGW-0100m_SBE39-5824_20140309.asc'
        sno{2} = '5824';
        plannedz{2} = 100;

        inst{3} = 'SBE37';
        fname{3} = 'SGW-0200m_SBE37-2952_20140311.asc'
        sno{3} = '2952';
        plannedz{3} = 200;

        inst{4} = 'SBE39';
        fname{4} = 'SGW-0400m_SBE39-5825_20140309.asc'
        sno{4} = '5825';
        plannedz{4} = 400;

        inst{5} = 'SBE37';
        fname{5} = 'SGW-0700m_SBE37-2954_20140311.asc'
        sno{5} = '2954';
        plannedz{5} = 700;



    case 2

                % ST Georges East
        moor_name = 'StGeorgesEast';
        mn = 'SGE'


        lat = -4-6.174/60;
        lon = 152+31.116/60;
        waterz = 1433;

        % deployment times
        start_time=datenum(2012,07,20,06,10,00);
        end_time = datenum(2014,03,8,04,00,00);

        % need to get rid of headers on SBE56 csv files except date,time,temp

        inst{1} = 'SBE37';
        fname{1} = 'SGE-0180m_SBE37-8715_20140311.asc'
        sno{1} = '8715';
        plannedz{1} = 180;


        inst{2} = 'SBE56';
        fname{2} = 'SGE-0200m_SBE56-00947_2014-03-09.asc.csv'
        sno{2} = '947';
        plannedz{2} = 200;

        inst{3} = 'SBE37';
        fname{3} = 'SGE-0300m_SBE37-8716_20140311.asc'
        sno{3} = '8716';
        plannedz{3} = 300;

        inst{4} = 'SBE39';
        fname{4} = 'SGE-0500m_SBE39-1472_20140309.asc'
        sno{4} = '1472';
        plannedz{4} = 500;

        inst{5} = 'SBE56';
        fname{5} = 'SGE-0600m_SBE56-00948_2014-03-09.asc.csv'
        sno{5} = '948';
        plannedz{5} = 600;
        
        inst{6} = 'SBE37';
        fname{6} = 'SGE-0700m_SBE37-8717_20140311.asc'
        sno{6} = '8717';
        plannedz{6} = 700;
        
    case 3
        
                % Vitiaz WEst
        moor_name = 'VitiazWest';
        mn = 'VW'


        lat = -5-58.69/60;
        lon = 147+39.97/60;
        waterz = 980;

        % deployment times
        start_time=datenum(2012,07,28,06,15,00);
        end_time = datenum(2014,03,13,18,45,00);

        % need to get rid of headers on SBE56 csv files except date,time,temp


        inst{1} = 'SBE39';
        fname{1} = 'VW-0250m_SBE039-1328_20140314.asc';
        sno{1} = '1328';
        plannedz{1} = 250;
        
        inst{2} = 'SBE39';
        fname{2} = 'VW-0400m_SBE039-1329_20140314.asc';
        sno{2} = '1329';
        plannedz{2} = 400;
        
        

        
        inst{3} = 'SBE56';
        fname{3} = 'VW-0500m_SBE56-00931_20140314.asc.csv'
        sno{3} = '931';
        plannedz{3} = 500;

        
        inst{4} = 'SBE56';
        fname{4} = 'VW-0650m_SBE56-00932_20140314.asc.csv'
        sno{4} = '932';
        plannedz{4} = 650;

        inst{5} = 'SBE37';
        fname{5} = 'VW-0750m_SBE37-8720_20140316.asc'
        sno{5} = '8720';
        plannedz{5} = 750;
        
        
    case 4
        moor_name = 'VitiazMiddle';
        mn = 'VM';
        
           lat = -5-56.65/60;
        lon = 147+46.68/60;
        waterz = 1131;
        
        start_time=datenum(2012,07,29,01,00,00);
        end_time = datenum(2014,03,14,02,00,00);

        % need to get rid of headers on SBE56 csv files except
        % date,time,temp

        inst{1} = 'SBE56';
        fname{1} = 'VM-0060m_SBE56-00933_20140314.asc.csv'
        sno{1} = '933';
        plannedz{1} = 60;
        

        inst{2} = 'SBE56';
        fname{2} = 'VM-0080m_SBE56-00934_20140314.asc.csv'
        sno{2} = '934';
        plannedz{2} = 80;     

                inst{3} = 'SBE56';
        fname{3} = 'VM-0130m_SBE56-00935_20140314.asc.csv'
        sno{3} = '935';
        plannedz{3} = 130;  
        
                inst{4} = 'SBE56';
        fname{4} = 'VM-0160m_SBE56-00936_20140314.asc.csv'
        sno{4} = '936';
        plannedz{4} = 160;  
        
                inst{5} = 'SBE56';
        fname{5} = 'VM-0250m_SBE56-00937_20140314.asc.csv'
        sno{5} = '937';
        plannedz{5} = 250;  
        
                inst{6} = 'SBE56';
        fname{6} = 'VM-0350m_SBE56-00938_20140314.asc.csv'
        sno{6} = '938';
        plannedz{6} = 350;  
        
                inst{7} = 'SBE56';
        fname{7} = 'VM-0400m_SBE56-00939_20140314.asc.csv'
        sno{7} = '939';
        plannedz{7} = 400;  
        
                inst{8} = 'SBE56';
        fname{8} = 'VM-0700m_SBE56-00928_20140314.asc.csv'
        sno{8} = '928';
        plannedz{8} = 700;  
        
                inst{9} = 'SBE56';
        fname{9} = 'VM-0650m_SBE56-00940_20140314.asc.csv'
        sno{9} = '940';
        plannedz{9} = 650;  

                
        inst{10} = 'SBE39';
        fname{10} = 'VM-0100m_SBE039-5820_20140314.asc'
        sno{10} = '5820';
        plannedz{10} = 100;  
        
                
        inst{11} = 'SBE39';
        fname{11} = 'VM-0180m_SBE039-5821_20140314.asc'
        sno{11} = '5821';
        plannedz{11} = 180;  
        
        inst{12} = 'SBE39';
        fname{12} = 'VM-0500m_SBE039-5822_20140314.asc'
        sno{12} = '5822';
        plannedz{12} = 500; 
        
        inst{13} = 'SBE37';
        fname{13} = 'VM-0150m_SBE37-8721_20140316.asc'
        sno{13} = '8721';
        plannedz{13} = 150; 

                
        inst{14} = 'SBE37';
        fname{14} = 'VM-0200m_SBE37-8722_20140316.asc'
        sno{14} = '8722';
        plannedz{14} = 200; 
        
                
        inst{15} = 'SBE37';
        fname{15} = 'VM-0300m_SBE37-8723_20140316.asc'
        sno{15} = '8723';
        plannedz{15} = 300; 
        
                
        inst{16} = 'SBE37';
        fname{16} = 'VM-0800m_SBE37-8724_20140316.asc'
        sno{16} = '8724';
        plannedz{16} = 800; 
        

    case 5                   % Vitiaz East
        moor_name = 'VitiazEast';
        mn = 'VE'


        lat = -5-54.98/60;
        lon = 147+50.5/60;
        waterz = 898;

        % deployment times

        
        start_time=datenum(2012,07,28,01,30,00);
        end_time = datenum(2014,03,14,21,00,00);

        % need to get rid of headers on SBE56 csv files except
        % date,time,temp

        inst{1} = 'SBE56';
        fname{1} = 'VE-0500m_SBE56-00944_20140315.asc.csv'
        sno{1} = '944';
        plannedz{1} = 500;

        
        inst{2} = 'SBE56';
        fname{2} = 'VE-0750m_SBE56-00946_20140315.asc.csv'
        sno{2} = '946';
        plannedz{2} = 750;

        
        inst{3} = 'SBE37';
        fname{3} = 'VE-0650m_SBE37-2951_20140315.asc'
        sno{3} = '2951';
        plannedz{3} = 650;
        
end

%% read data


nf = length(fname);

for inf = 1:nf

    clear t
    infile = [main_path file_dat moor_name '/' char(inst{inf}) '/' fname{inf}]


    % SBE56 temp only
    if (inst{inf}=='SBE56')

        display('reading SBE56')
        [time,temp] = importsBE56(infile);
        display('finished reading SBE56')

        t.time = time;
        t.temperature = temp;
        t.pressure = zeros(length(temp),1);

    end



    % default pressure and read pass back conductivity
  if (inst{inf}=='SBE37' | inst{inf}=='SBE39')
        t = readSBE37(infile,lat,lon,plannedz{inf},1);
    end

    % save data within start and end times
    isub = (t.time>=start_time&t.time<=end_time);
    t.time = t.time(isub);
    t.temperature = t.temperature(isub);
    t.pressure = t.pressure(isub);

    if (inst{inf}=='SBE37')
        t.salinity = t.salinity(isub);
        t.conductivity = t.conductivity(isub);
    end

    t.startdate_mat = t.time(1);
    t.stopdate_mat = t.time(end);
    t.latitude = lat;
    t.longitude = lon;
    t.plannedMeterDepth = plannedz{inf};
    t.waterDepth = waterz;
    t.meterType = inst{inf};
    t.serialNo = sno{inf};


%% plot raw data
    figure(1)
    clf

    set(gcf,'defaultaxesFontWeight','bold')
    set(gcf,'defaultAxesFontSize',15)
    set(gcf,'DefaultAxesLineWidth',2)

    subplot(3,1,1)
    plot(t.time,t.temperature)
    axis([t.time(1) t.time(end) min(t.temperature)-1 max(t.temperature)+1])
    %datetick('x','m','keeplimits');
    montick('x','m',get(gca,'position'))
    ylabel('Temperature (deg C)')
    title([t.meterType ' ' t.serialNo ' ' moor_name ' z= ' num2str(t.plannedMeterDepth)])

    if (t.meterType=='SBE37')
        subplot(3,1,3)
        plot(t.time,t.salinity)
        axis([t.time(1) t.time(end) min(t.salinity)-1 max(t.salinity)+1])
        axis([t.time(1) t.time(end) 33 35])
        %datetick('x',12,'keeplimits');
        montick('x','m',get(gca,'position'))

        ylabel('Salinity (psu)')
    end


    subplot(3,1,2)
    plot(t.time,t.pressure)
    axis([t.time(1) t.time(end) min(t.pressure)-5 max(t.pressure)+5])
    %datetick('x',12,'keeplimits');
    ylabel('Pressure (db)')
    montick('x','m',get(gca,'position'))
    outfile = [main_path file_fig moor_name '/' moor_name '_' inst{inf} '_' sno{inf} '.pdf']
    eval(['print -dpdf ' outfile]);
    %% save data

    tl_file = [main_path file_dat moor_name '/' inst{inf} '/' moor_name '_' inst{inf} '_' sno{inf} '_qc.mat'];

    eval(['save ' tl_file ' t']);
    

    %%
end % next data file
