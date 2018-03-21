% adcp_qc.m

% Programmer: Janet Sprintall
% Date: 5 October 2013

% Purpose: program to take the _qc.mat adcp data sets from each mooring and
% a) filter the data to the required time interval
% b) interpolate the data to a common time interval
% c) assign pressure to those instruments on the mooring without pressure
% d) save the data for each mooring as a stacked interval

% had thought about including aquadopp temperature but only 3 m from most
% shallowest and too tricky to do so discounted.

close all
clear all

%addpath(path,'/Applications/MATLAB74/toolbox/plots')

% common paths

main_path = '/Users/jsprintall/Documents/JanetsDocuments/southernocean/AntarcticDipoleMode/Mooring/';
% moor_deployment


%%

start_time=datenum(2011,04,11,23,00,00);
end_time = datenum(2013,03,24,15,30,00);    % have last good data at 27-Jun-2007 21:30:35
% for file names
moor_name = 'ADP_';

% path for prints and other I/O in process subdirectory
moor_path = 'Data/SBE/';


% strip out the _ for the figure titles as it becomes subscript
mn = strrep(moor_name,'_',' ');

samplerate = 30;    % choose sample rate for filtering purposes

% Interpolate to common time basis
% using start_time and end_time
% first set up a common time axis
% using function dateax
% starting at start_time, with 1/2 hour sampling interval (1\48 hours) and
% with an array length from the endtime-starttime in hours

ndate = dateax(datevec(start_time),1/48,48*(end_time-start_time)+1);

%% 2. READ IN THE DATA FILES and sort mooring from shallowest to deepest

% load in all adcp qc data for that mooring and deployment
mdir = dir([main_path moor_path moor_name '*_qc.mat']);

% sort out order via planned meter depth
nfiles =0;

for inst = 1:length(mdir);
    clear t
    eval(['load ', main_path,moor_path,mdir(inst).name]);
    if exist('t')
        nfiles = nfiles+1;
        mdepth(nfiles) = (t.plannedMeterDepth);
        eval(['t' t.serialNo '= t;']);
        % keep a cell structure of the instrument files and other goodies
        isn{nfiles} = t.serialNo;
        m_depth{nfiles} = mdepth(nfiles);
        mtype{nfiles} = t.meterType;
        t.serialNo
        length(t.time)
    end
end

% okay, now just sort out those instruments that exist and put them in
% order from shallowest to deepest;
[junk,iz] = sort(mdepth);
mdir = mdir(iz);
isn = isn(iz);
mdepth = mdepth(iz);
mtype = mtype(iz);

% remember that deepest instrment is at ~4600 m and is a SBE39P so does not
% need pressure and is probably not good to include here

%% stack data sampleinterval time step   and fill in missing depths for
% temperature

got1 = 0;

for is = 1:nfiles
    clear t isub
    eval(['t = t' char(isn(is)) ';']);
    t.serialNo
    isub = find(isnan(t.pressure));
    if (length(isub) ~= length(t.pressure))
        % case where have pressure data so sort into ascending pressure
        % (positive down)
        % save the pressure and target pressure as this may be important
        display('got1 with pressure')
        got1 = is;
    else
        % no pressure data so use the data from the shallower instrument
        if (got1)
            % read in pressure data from above instrument

            eval(['t_before = t' char(isn(got1)) ';']);
            readme_p = ['pressure record ' (t.serialNo) ' from ' (t_before.serialNo)];

            % this is the planned wire distance between these instruments

            zdiff = t.plannedMeterDepth - t_before.plannedMeterDepth;
% this is the interpolated pressure for the missing data 
            
            pressure = t_before.pressure + zdiff;
            t.pressure = pressure;
%            eval(['t.pressure = t' char(isn(got1)) ';']);
                eval(['t' char(isn(is)) '= t;']);

            
        end
    end

end


%% save data in structure

[T, S, pressure] = deal(nan(length(ndate),nfiles));

% for pressure and temperature need 2 fields
for ii =1:nfiles
    clear t
    eval(['t = t' char(isn(ii)) ';']);
    T(:,ii) = t.temperature;
    pressure(:,ii) = t.pressure;
    
    % now if salinity save that too
    if(exist('t.salinity'))
        S(:,ii) = t.salinity;
    end
    
end

lat = t.latitude;
lon = t.longitude;
waterz = t.waterDepth;
%moor_name = s.mooringName;

%%
t =struct('temperature',T);
t =setfield(t,'pressure',pressure);
t =setfield(t,'salinity',S);
t =setfield(t,'time',ndate);
t =setfield(t,'latitude',lat);
t =setfield(t,'longitude',lon);
t =setfield(t,'meterType',mtype);
t =setfield(t,'serialNo',isn);
t =setfield(t,'waterDepth',waterz);
t =setfield(t,'plannedMeterDepth',m_depth);
t =setfield(t,'sampleRate',samplerate);
t =setfield(t,'mooringName',moor_name);


display('saving file')
T_file = [main_path moor_path moor_name '_allSBE.mat'];
eval(['save ' T_file ' t']);

%% 6. PLOT THE temperature DATA

display('plotting data')

[nt,nz] = size(t.temperature);

figure(103)
clf
plot(ndate,t.temperature(:,1:nz-1))
shading flat
colormap('bluewhitered')
        montick('x','m',get(gca,'position'))
ylabel('Temperature')
eval(['print -djpeg ',[main_path 'Figures/' moor_name],'_allSBE.jpg'])
