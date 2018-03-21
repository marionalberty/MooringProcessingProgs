% NB: read data (*.dat file) using Import data function
% this program just reads in data and saves as structure

clear all
close all
%addpath('/Applications/MATLAB74/toolbox/plots')
% main_path = '/Users/jsprintall/Documents/JanetsDocuments/swpacific/MoorSPICE/Moorings/';
main_path = '/Users/marionsofiaalberty/MATLAB/Solomon_Sea/Moorings/';

moor_depl = 8;

switch moor_depl
  case 1
    % 1. St Georges East
    % --------------
    % ---------------
    start_time=datenum(2012,07,20,06,10,00);
    end_time = datenum(2014,03,8,04,00,00);
    % for file names
    moor_name = 'StGeorgesEast';
    
    lat = -4-6.174/60;
    lon = 152+31.116/60;
    waterz = 1433;
    
    % path for prints and other I/O in process subdirectory
    moor_path = ['Data/Raw/' moor_name '/AquaDopp/'];
    ShortName = 'SGE';
    NoAQD = 1;
    deployDepth{1} = 750;
    SerNo{1} = 9906;
    MagDec = 0.0;
    % needs to be a cell??            DownloadDate(1) = '20140309_01';
    
    MeasurementInterval = 3600; %(sec)
    MeasurementInterval_unit = 'seconds';
    
  case 2
    % Vitiaz West
    moor_name = 'VitiazWest';
    % path for prints and other I/O in process subdirectory
    moor_path = ['Data/Raw/' moor_name '/AquaDopp/'];
    ShortName = 'VW';
    NoAQD = 1;
    deployDepth{1} = 800;
    SerNo{1} = 9908;
    % needs to be a cell??            DownloadDate(1) = '20140309_01';
    
    waterz = 980;
    MagDec = 0.0;
    
    MeasurementInterval = 3600; %(sec)
    MeasurementInterval_unit = 'seconds';
    
    lat = -5-58.69/60;
    lon = 147+39.97/60;
    
    % deployment times
    start_time=datenum(2012,07,28,06,30,00);
    end_time = datenum(2014,03,13,20,00,00);
    
  case 3
    % Vitiaz Middle
    moor_name = 'VitiazMiddle';
    % path for prints and other I/O in process subdirectory
    moor_path = ['Data/Raw/' moor_name '/AquaDopp/'];
    ShortName = 'VM';
    NoAQD = 2;
    deployDepth{1} = 850;
    SerNo{1} = 6832;
    
    deployDepth{2} = 1000;
    SerNo{2} = 9902;
    % needs to be a cell??            DownloadDate(1) = '20140309_01';
    
    waterz = 1130;
    MagDec = 0.0;
    
    MeasurementInterval = 3600; %(sec)
    MeasurementInterval_unit = 'seconds';
    
    lat = -5-56.65/60;
    lon = 147+46.68/60;
    
    % deployment times
    start_time=datenum(2012,07,29,01,00,00);
    end_time = datenum(2014,03,14,02,00,00);
    
    
  case 4
    % Vitiaz East
    moor_name = 'VitiazEast';
    
    % path for prints and other I/O in process subdirectory
    moor_path = ['Data/Raw/' moor_name '/AquaDopp/'];
    ShortName = 'VE';
    NoAQD = 1;
    deployDepth{1} = 700;
    SerNo{1} = 9931;
    
    waterz = 900;
    MagDec = 0.0;
    
    MeasurementInterval = 3600; %(sec)
    MeasurementInterval_unit = 'seconds';
    
    lat = -5-54.98/60;
    lon = 147+50.05/60;
    
    % deployment times
    start_time=datenum(2012,07,28,01,30,00);
    end_time = datenum(2014,03,15,20,00,00);
    
  case 5
    % Solomon M1
    moor_name = 'Solomon_M1';
    % path for prints and other I/O in process subdirectory
    moor_path = ['Data/Raw/' moor_name '/AquaDopp/'];
    ShortName = 'M1';
    NoAQD = 2;
    deployDepth{1} = 100;
    SerNo{1} = 5939;
    deployDepth{2} = 500;
    SerNo{2} = 1959;
    
    waterz = 2050;
    MagDec = 6.39;
    
    MeasurementInterval = 1800; %(sec)
    MeasurementInterval_unit = 'seconds';
    
    lat = -4-57.408/60;
    lon = 153+16.149/60;
    
    % deployment times
    start_time=datenum(2012,07,22,00,00,00);
    end_time = datenum(2014,03,07,06,20,03);
    
  case 6
    % Solomon M2a
    moor_name = 'Solomon_M2a';
    % path for prints and other I/O in process subdirectory
    moor_path = ['Data/Raw/' moor_name '/AquaDopp/'];
    ShortName = 'M2a';
    NoAQD = 1;
    deployDepth{1} = 200;
    SerNo{1} = 2295;
    
    waterz = 2560;
    MagDec = 6.39;
    
    MeasurementInterval = 1800; %(sec)
    MeasurementInterval_unit = 'seconds';
    
    lat = -5-09.424/60;
    lon = 153+17.016/60;
    
    % deployment times
    start_time=datenum(2012,07,22,00,00,00);
    end_time = datenum(2014,03,07,06,00,00);
    
  case 7
    % Solomon M2b
    moor_name = 'Solomon_M2b';
    % path for prints and other I/O in process subdirectory
    moor_path = ['Data/Raw/' moor_name '/AQUADOPP/'];
    ShortName = 'M2b';
    NoAQD = 1;
    deployDepth{1} = 500;
    SerNo{1} = 2323;
    
    waterz = 2710;
    MagDec = 6.39;
    
    MeasurementInterval = 1800; %(sec)
    MeasurementInterval_unit = 'seconds';
    
    lat = -5-09.312/60;
    lon = 153+19.686/60;
    
    % deployment times
    start_time=datenum(2012,07,22,00,00,00);
    end_time = datenum(2014,03,07,05,30,00);
    
  case 8
    % Solomon M3
    moor_name = 'Solomon_M3';
    % path for prints and other I/O in process subdirectory
    moor_path = ['Data/Raw/' moor_name '/AQUADOPP/'];
    ShortName = 'M3';
    NoAQD = 2;
    deployDepth{1} = 100;
    SerNo{1} = 1964;
    deployDepth{2} = 500;
    SerNo{2} = 2325;
    
    waterz = 2617;
    MagDec = 6.39;
    
    MeasurementInterval = 3600; %(sec)
    MeasurementInterval_unit = 'seconds';
    
    lat = -5-08.447/60;
    lon = 154+17.81/60;
    
    % deployment times
    start_time=datenum(2012,07,22,00,00,00);
    end_time = datenum(2014,03,05,03,57,04);
end

for AQD = 1:NoAQD;
  
  fname = [main_path moor_path ShortName '_' num2str(deployDepth{AQD})...
    'm_AQD_' num2str(SerNo{AQD}) '.mat'];
  eval(['load ' fname]);
  % eval(['AQDraw = ' ShortName '_' num2str(deployDepth{AQD}) 'm_AQD_'...
  % num2str(SerNo{AQD}) ';']);
  
  
  
  %% now set variables
  % %
  % % 1   Month                            (1-12)
  % %  2   Day                              (1-31)
  % %  3   Year
  % %  4   Hour                             (0-23)
  % %  5   Minute                           (0-59)
  % %  6   Second                           (0-59)
  % %  7   Error code
  % %  8   Status code
  % %  9   Velocity (Beam1|X|East)          (m/s)
  % % 10   Velocity (Beam2|Y|North)         (m/s)
  % % 11   Velocity (Beam3|Z|Up)            (m/s)
  % % 12   Amplitude (Beam1)                (counts)
  % % 13   Amplitude (Beam2)                (counts)
  % % 14   Amplitude (Beam3)                (counts)
  % % 15   Battery voltage                  (V)
  % % 16   Soundspeed                       (m/s)
  % % 17   Soundspeed used                  (m/s)
  % % 18   Heading                          (degrees)
  % % 19   Pitch                            (degrees)
  % % 20   Roll                             (degrees)
  % % 21   Pressure                         (dbar)
  % % 22   Pressure                         (m)
  % % 23   Temperature                      (degrees C)
  % % 24   Analog input 1
  % % 25   Analog input 2
  % % 26   Speed                            (m/s)
  % % 27   Direction                        (degrees)
  % -----------------------------------------------------
  
  t1 = 50;
  t2 = length(AQDraw)-56;
  
  Month = AQDraw(t1:t2,1);
  Day  = AQDraw(t1:t2,2);
  Year= AQDraw(t1:t2,3);
  Hour    = AQDraw(t1:t2,4);
  Minute= AQDraw(t1:t2,5);
  Second   = AQDraw(t1:t2,6);
  
  % convert to matlab time
  time = datenum(Year,Month,Day,Hour,Minute,Second);
  
  
  %  7   Error code
  ErrorCode = AQDraw(t1:t2,7);
  %  8   Status code
  StatusCode = AQDraw(t1:t2,7);
  
  UVelocity = AQDraw(t1:t2,9);  %     (m/s)
  VVelocity = AQDraw(t1:t2,10);%        (m/s)
  WVelocity = AQDraw(t1:t2,11);  %     (m/s)
  
  UVelocity_unit = 'm/s';
  VVelocity_unit = 'm/s';
  WVelocity_unit = 'm/s';
  
  
  
  % % 12   Amplitude (Beam1)                (counts)
  % % 13   Amplitude (Beam2)                (counts)
  % % 14   Amplitude (Beam3)                (counts)
  % %
  % % 15   Battery voltage                  (V)
  % % 16   Soundspeed                       (m/s)
  AmplitudeBeam1 = AQDraw(t1:t2,12);
  AmplitudeBeam2 = AQDraw(t1:t2,13);
  AmplitudeBeam3 = AQDraw(t1:t2,14);
  AmplitudeUnit = 'counts';
  
  BatteryVoltage = AQDraw(t1:t2,15);
  BatteryVoltage_units = 'V';
  
  Soundspeed = AQDraw(t1:t2,17); %used                  (m/s)
  Soundspeed_unit = 'm/s';
  
  Heading = AQDraw(t1:t2,18);%                          (degrees)
  Pitch = AQDraw(t1:t2,19);%                           (degrees)
  Roll = AQDraw(t1:t2,20);%                            (degrees)
  
  Heading_unit = 'degrees';
  Pitch_unit = 'degrees';
  Roll_unit = 'degrees';
  
  PressureDb = AQDraw(t1:t2,21);%                        (dbar)
  PressureM = AQDraw(t1:t2,22);%                    (m)
  
  PressureDb_unit = 'dbar';
  PressureM_unit = 'm';
  
  Temperature = AQDraw(t1:t2,23);%                     (degrees C)
  
  Temperature_unit = 'degrees C';
  Speed = AQDraw(t1:t2,26);%                            (m/s)
  Direction = AQDraw(t1:t2,27);%                        (degrees)
  
  Speed_unit = 'm/s';
  Direction_unit = 'degrees';
  
  %% now save only that part that was during the mooring deployment
  
  isub = find(time>=datenum(start_time)&time<=datenum(end_time));
  time = time(isub);
  UVelocity = UVelocity(isub);
  VVelocity = VVelocity(isub);
  WVelocity = WVelocity(isub);
  Speed = Speed(isub);
  Direction = Direction(isub);
  Soundspeed = Soundspeed(isub);
  Heading = Heading(isub);
  Pitch = Pitch(isub);
  Roll = Roll(isub);
  PressureDb = PressureDb(isub);
  Pressurem = PressureM(isub);
  Temperature = Temperature(isub);
  
  AmplitudeBeam1 = AmplitudeBeam1(isub);
  AmplitudeBeam2 = AmplitudeBeam2(isub);
  AmplitudeBeam3 = AmplitudeBeam3(isub);
  
  BatteryVoltage = BatteryVoltage(isub);
  
  %% now save to outfile
  t =struct('time',time);
  
  t = setfield(t,'mooring',moor_name);
  
  t = setfield(t,'latitude', lat);
  t = setfield(t,'longitude', lon);
  t = setfield(t,'waterDepth', waterz);
  
  t =setfield(t,'sampleRate',MeasurementInterval);
  t =setfield(t,'sampleRate_unit',MeasurementInterval_unit);
  
  t = setfield(t,'u',UVelocity');
  t = setfield(t,'v',VVelocity');
  t = setfield(t,'w',WVelocity');
  t = setfield(t,'speed',Speed');
  t = setfield(t,'direction',Direction');
  
  
  t = setfield(t,'velocity_unit',UVelocity_unit);
  t = setfield(t,'direction_unit','degrees');
  
  t = setfield(t,'soundspeed',Soundspeed);
  t = setfield(t,'soundspeed_unit',Soundspeed_unit);
  
  t = setfield(t,'heading',Heading');
  t = setfield(t,'pitch',Pitch');
  t = setfield(t,'roll',Roll');
  
  t = setfield(t,'amplitudeBeam1',AmplitudeBeam1');
  t = setfield(t,'amplitudeBeam2',AmplitudeBeam2');
  t = setfield(t,'amplitudeBeam3',AmplitudeBeam3');
  t = setfield(t,'amplitudeUnit',AmplitudeUnit);
  t = setfield(t,'batteryVoltage',BatteryVoltage');
  t = setfield(t,'batteryVoltageUnits',BatteryVoltage_units);
  
  t = setfield(t,'heading_unit',Heading_unit);
  t = setfield(t,'pitch_unit',Pitch_unit);
  t = setfield(t,'roll_unit',Roll_unit);
  
  t = setfield(t,'pressure_dbar',PressureDb');
  t = setfield(t,'pressure_meters',PressureM');
  
  
  t = setfield(t,'temperature',Temperature');
  
  
  t = setfield(t,'targetDeploymentDepth',deployDepth{AQD});
  t = setfield(t,'serialNumber',SerNo{AQD});
  
  t = setfield(t,'magneticDeclination',MagDec);
  t = setfield(t,'magDec_README',...
    'Velocity and Speed are UNCORRECTED for magdec (degE)');
  
  
  fname = [main_path moor_path moor_name '_' ...
    num2str(t.targetDeploymentDepth) '_' num2str(t.serialNumber) ...
    '_AQD_qc.mat'];
  eval(['save ' fname ' t'])
  
  %% now plot some of the key variables
  
  %     display('plotting data')
  %
  %
  %     figure(103)
  %     clf
  %     subplot(2,1,1)
  %     plot(t.time,t.U)
  %     title([ShortName ' ' num2str(t.SerialNumber) ...
  %     ' AquaDopp Velocity (Mean Pressure) ' num2str(nanmean(t.Pressure_dbar))])
  %     ylabel('U (m/s)')
  %     hold on
  %     plot([t.time(1) t.time(end)],[0 0],'r-')
  %     axis([t.time(1) t.time(end) -inf inf])
  %     montick('x','m',get(gca,'position'))
  %
  %     subplot(2,1,2)
  %     plot(t.time,t.V)
  %     ylabel('V (m/s)')
  %     hold on
  %     plot([t.time(1) t.time(end)],[0 0],'r-')
  %     axis([t.time(1) t.time(end) -inf inf])
  %     montick('x','m',get(gca,'position'))
  %
  %     eval(['print -dpdf ',[main_path 'Figures/' moor_name '/'   ...
  %     moor_name '_' num2str(t.targetDeploymentDepth) 'm_' ...
  %     num2str(t.SerialNumber)],'_AQD_Velocity.pdf'])
  %
  %
  %     figure(13)
  %     clf
  %     subplot(3,1,1)
  %     plot(t.time,t.speed)
  %     title([ShortName ' ' num2str(t.SerialNumber) ' AquaDopp Velocity (Mean Pressure) ' num2str(nanmean(t.Pressure_dbar))])
  %     ylabel('Speed (m/s)')
  %     axis([t.time(1) t.time(end) -inf inf])
  %     montick('x','m',get(gca,'position'))
  %
  %     subplot(3,1,2)
  %     plot(t.time,t.direction,'.')
  %     ylabel('Direction (degrees)')
  %     axis([t.time(1) t.time(end) -inf inf])
  %     montick('x','m',get(gca,'position'))
  %
  %
  %     eval(['print -dpdf ',[main_path 'Figures/' moor_name '/'   moor_name '_' num2str(t.targetDeploymentDepth) 'm_' num2str(t.SerialNumber)],'_AQD_speed.pdf'])
  %
  %     figure(1)
  %     clf
  %     subplot(3,1,1)
  %     plot(t.time,t.Pressure_dbar)
  %     title([ShortName ' ' num2str(t.SerialNumber) ' AquaDopp Velocity (Mean Pressure) ' num2str(nanmean(t.Pressure_dbar))])
  %     ylabel('Pressure (db)')
  %     axis([t.time(1) t.time(end) -inf inf])
  %     montick('x','m',get(gca,'position'))
  %
  %     subplot(3,1,2)
  %     plot(t.time,t.Temperature)
  %     ylabel('Temperature (degrees)')
  %     axis([t.time(1) t.time(end) -inf inf])
  %     montick('x','m',get(gca,'position'))
  %
  %
  %     eval(['print -dpdf ',[main_path 'Figures/' moor_name '/'   moor_name '_' num2str(t.targetDeploymentDepth) 'm_' num2str(t.SerialNumber)],'_AQD_Velocity_Pressure_T.pdf'])
  %
  %
  %
  %     figure(10)
  %     clf
  %     subplot(4,1,1)
  %     plot(t.time,t.Heading,'.')
  %     title([ShortName ' ' num2str(t.SerialNumber) ' AquaDopp Diagnostics (Mean Pressure) ' num2str(nanmean(t.Pressure_dbar))])
  %     ylabel('heading (deg)')
  %     axis([t.time(1) t.time(end) -inf inf])
  %     montick('x','m',get(gca,'position'))
  %
  %     subplot(4,1,2)
  %     plot(t.time,t.Pitch,'.')
  %     ylabel('pitch (deg)')
  %     axis([t.time(1) t.time(end) -inf inf])
  %     montick('x','m',get(gca,'position'))
  %
  %     subplot(4,1,3)
  %     plot(t.time,t.Roll,'.')
  %     ylabel('roll (deg)')
  %     axis([t.time(1) t.time(end) -inf inf])
  %     montick('x','m',get(gca,'position'));
  %
  %     subplot(4,1,4)
  %     plot(t.time,t.BatteryVoltage,'.')
  %     ylabel('Battery Voltage (V)')
  %     axis([t.time(1) t.time(end) -inf inf])
  %     montick('x','m',get(gca,'position'));
  %
  %     eval(['print -dpdf ',[main_path 'Figures/' moor_name '/'   moor_name '_' num2str(t.targetDeploymentDepth) 'm_' num2str(t.SerialNumber)],'_AQD_Diagnostics.pdf'])
  
end % done for each aquadopp