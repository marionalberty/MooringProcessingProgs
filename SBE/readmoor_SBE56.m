% read in SBEs from mooring data 

clear all
close all

addpath('/Applications/MATLAB74/toolbox/plots')


main_path = '/Users/jsprintall/Documents/JanetsDocuments/swpacific/MoorSPICE/Moorings/';

moor_depl = 1;

%%
% 1. INITIALIZE INPUT FOR EACH MOORING
%
% These parameters probably need to be changed for each mooring and
% deployment!

switch moor_depl
 case 1
     % ST Geirges West
     moor_name = 'StGeorgesWest';
     mn = 'SGE'
% path for prints and other output in process subdirectory
filen = 'Figures/'
% path for input file from raw subdirectory

lat = -4-6.82/60;
lon = 152+31.116/60;
waterz = 1243;

% deployment times
start_time=datenum(2012,07,19,23,30,00);
end_time = datenum(2014,03,14,01,00,00);

% cell arrays of SBE37


inst{1} = 'SBE56';
fname{1} = 'SGW-0600m_SBE56-01024_2014-03-09.asc.csv'
sno{1} = '01024';
plannedz{1} = 600;

% % inst{1} = 'SBE37';
% % fname{1} = 'sbe37-2949.asc'
% % sno{1} = '2949';
% % plannedz{1} = 250;
% % 
% % inst{2} = 'SBE37';
% % fname{2} = 'sbe37-2951.asc'
% % sno{2} = '2951';
% % plannedz{2} = 350;
% % 
% % inst{3} = 'SBE39';
% % fname{3} = 'sbe39-1326.asc'
% % sno{3} = '1326';
% % plannedz{3} = 300;
% % 
% % inst{4} = 'SBE39';
% % fname{4} = 'sbe39-1327.asc'
% % sno{4} = '1327';
% % plannedz{4} = 450;

case 2
     % North Mindoro
     moor_name = 'Mindoro_deploy2_';
% path for prints and other output in process subdirectory
filen = 'deploy2/nmin/process/'
% path for input file from raw subdirectory
filenr = 'deploy2/nmin/raw/'

lat = 11.8933;
lon = 121.0549;
waterz = 450;

% deployment times
start_time=datenum(2007,12,21,23,30,00);
% at the moment save out past the retrieval to determine what the mooring
% did
%end_time = datenum(2009,03,19,11,00,00);

% now real end of good data time
end_time = datenum(2009,3,18,22,00,00);

% cell arrays of SBE37

inst{1} = 'SBE37';
fname{1} = 'sbe37-2952.asc'
sno{1} = '2952';
plannedz{1} = 200;

inst{2} = 'SBE37';
fname{2} = 'sbe37-2954.asc'
sno{2} = '2954';
plannedz{2} = 300;

inst{3} = 'SBE39';
fname{3} = 'sbe39-1328.asc'
sno{3} = '1328';
plannedz{3} = 250;

inst{4} = 'SBE39';
fname{4} = 'sbe39-1329.asc'
sno{4} = '1329';
plannedz{4} = 350;

inst{5} = 'SBE39';
fname{5} = 'sbe39-1471.asc'
sno{5} = '1471';
plannedz{5} = 443;

    case 3
        % surigao mooring deploy 2
 moor_name = 'Surigao_deploy2_';
% path for prints and other output in process subdirectory
filen = 'deploy2/suri/process/'
% path for input file from raw subdirectory
filenr = 'deploy2/suri/raw/'

lat = 10+26.063/60;
lon = 125+22.543/60;

waterz = 630;

% deployment times
start_time=datenum(2007,12,24,23,30,00);
% now real end of good data time
end_time = datenum(2009,3,18,22,00,00);

% cell arrays of SBE37

inst{1} = 'SBE39';
fname{1} = 'sbe39-1472.asc'
sno{1} = '1472';
plannedz{1} = 550;

% not good data yet
% inst{2} = 'SBE39';
% fname{2} = 'sbe39-1004.asc'
% sno{2} = '1004';
% plannedz{2} = 620;

end

%% read data


nf = length(fname);

for inf = 1:nf
    
    clear t
infile = [main_path 'Data/' moor_name '/' char(inst{inf}) '/' fname{inf}]


% default pressure and read pass back conductivity

if (inst{inf}=='SBE56')
    nheaders = 13;

[date,time,temp] = importsBE56(infile);

%    t = importsBE56(infile,lat,lon,plannedz{inf},1);

shite

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

set(gcf,'defaultaxesFontWeight','bold')
set(gcf,'defaultAxesFontSize',15)
set(gcf,'DefaultAxesLineWidth',2)

subplot(3,1,1)
plot(t.time,t.temperature)
axis([t.time(1) t.time(end) min(t.temperature)-1 max(t.temperature)+1])
datetick('x',12,'keeplimits');
ylabel('Temperature (deg C)')
title([t.meterType ' ' t.serialNo ' ' moor_name ' z= ' num2str(t.plannedMeterDepth)])

if (t.meterType=='SBE37')
subplot(3,1,3)
plot(t.time,t.salinity)
axis([t.time(1) t.time(end) min(t.salinity)-1 max(t.salinity)+1])
datetick('x',12,'keeplimits');
ylabel('Salinity (psu)')
end


subplot(3,1,2)
plot(t.time,t.pressure)
axis([t.time(1) t.time(end) min(t.pressure)-5 max(t.pressure)+5])
datetick('x',12,'keeplimits');
ylabel('Pressure (db)')

outfile = [main_path filen moor_name sno{inf} 'SBE.jpg']
eval(['print -djpeg ' outfile]);
%% save data 

tl_file = [main_path filen moor_name sno{inf} '_qc.mat'];

eval(['save ' tl_file ' t']);

%%
end % next data file
