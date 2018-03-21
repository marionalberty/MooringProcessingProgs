% calibrate_L1

% Given sensor drifts in temperature and salinity, correct the available
% data and save the file as Level_1, for the Solomon Sea moorings

clear all; close all; clc

set(0,'defaultaxesfontsize',12,'defaultaxeslinewidth',0.7,...
  'defaultlinelinewidth',1,'defaultpatchlinewidth',0.7,...
  'defaultFigureColor','white')

%% Set script options

params.print=1; % 1= save figures
params.save=1; % 1= save qcLevel 1 data into 'mat' file
params.debug=0; % 1= debug mode
params.qcLevel=1; % define qcLevel type

%% Set paths and directories

% Personnal paths:
% driveName='/Users/cyrilgermineaud/Documents/MATLAB/';
% addpath([driveName 'Routines_Cyril'])
driveName='/Users/marionsofiaalberty/MATLAB/Solomon_Sea/';


%% Use user input to get calibration data
keeprunning='y';

% Double check path
prompt=['Is this the correct personal path(y/n): ' driveName '\n'];
driveYN=input(prompt,'s');
if strcmp(driveYN,'n')
  prompt='Please type the correct personal path \n';
  driveName=input(prompt,'s');
end

while strcmp(keeprunning,'y')
  % Get mooring name
  prompt=['Which mooring? Use the corresponding data folder name,' ...
    'e.g.Solomon_M1 or VitiazEast \n'];
  mooring=input(prompt,'s');
  
  % Get instrument type
  prompt='What instrument? e.g. SBE37, AQD, or RBR1050?\n';
  inst=input(prompt,'s');
  
  % Get serial number
  prompt='What is the serial number? Include leading zeros.\n';
  sernum=input(prompt,'s');
  
  % Get deployment depth
  prompt='What is the deployment depth? Including leading zeros.\n';
  depdep=input(prompt,'s');
  
  % Set input/output paths
  params.inpath =[driveName 'Moorings/Data/Level_0/' mooring '/'];
  params.outpath =[driveName 'Moorings/Data/Level_1/' mooring '/'];
  params.pathfig = [driveName 'Moorings/Figures/Level_1/' mooring '/'];
  params.dirpng = [params.pathfig 'PNG'];
  params.dirpdf = [params.pathfig 'PDF'];
  params.moorName = mooring; % define mooring name
  
  % Load file
  disp('Loading file')
  load([params.inpath strrep(mooring,'_','') '_' inst '_' sernum '_' depdep '_0.mat'])
  
  %% Temperature correction
  disp('Alright, let''s start with the temperature sensor...')
  
  % Get pre-deployment calibration date
  prompt = 'What is the pre-deployment calibration date? [YYYY MM DD]\n';
  dnum_preT = datenum(input(prompt));
  
  % Get post-deployment calibration date
  prompt = 'What is the post-deployment calibration date? [YYYY MM DD]\n';
  dnum_posT = datenum(input(prompt));
  
  % Offset = b * (residual / n)
  % b = number of days between pre-cruise calibration and the cast
  % n = number of days between pre- and post-cruise calibrations
  % residual = residual from calibration sheet
  temp0 = t.temperature;
  b = t.time - dnum_preT;
  n = dnum_posT - dnum_preT;
  
  % See if they have a drift or residual
  prompt = 'Do you have a drift or residual? (d/r)?\n';
  DR = input(prompt,'s');
  
  if strcmp(DR,'d')
    % Get temperature drift
    prompt = 'What is the temperature drift [deg C/year]?\n';
    residual = -n*input(prompt)/365;
  else
    % Get temperature residual
    prompt = 'What is the temperature residual [mdeg C] from calibration?\n';
    residual = input(prompt)/1000;
  end
  
  % Correct temperature drift
  offset = b * (residual / n);
  t.temperature = temp0 + offset;
  
  % Save calibration information
  t.Temp_calibration.residual = residual;
  t.Temp_calibration.pre_cal_date = datestr(dnum_preT);
  t.Temp_calibration.pos_cal_date = datestr(dnum_posT);
  
  % Plot temperature calibration
  if params.print == 1
    fig=figure(1);
    set(gcf,'PaperUnits','centimeters')
    xSize = 30; ySize = 21;
    xLeft = (30-xSize)/2; yTop = (21-ySize)/2;
    set(gcf,'Position',[xLeft yTop xSize*40 ySize*40])
    
    plot(t.time,t.temperature,'k',t.time,temp0,'b')
    legend('Calibrated','Level 0')
    title(sprintf('%s_%s_%s_%s_%d',mooring,inst,sernum,depdep,...
      params.qcLevel),'interpreter','None','fontsize',20)
    datetick('x','mmm-yy','keeplimits')
    ylabel('Temperature (°C)')
    grid on;box on
    
    set(gcf,'units','centimeters')
    set(gcf,'papersize',[32 20])
    set(gcf,'paperposition',[0,0,32,20])
    nfop=sprintf('%s_%s_%s_%s_%d',mooring,inst,sernum,depdep,...
      params.qcLevel);
    nfoppng=sprintf('%s/%s.png',params.dirpng,nfop);
    print(gcf,'-dpng','-r200',nfoppng);
    
    nfoppdf=sprintf('%s/%s.pdf',params.dirpdf,nfop);
    print(gcf,'-dpdf','-r200',nfoppdf);
    close
  end
  
  %% Conductivity correction
  
  prompt = 'Do you have conductivity calibration data? (y/n)\n';
  conYN = input(prompt,'s');
  if strcmp(conYN,'y')
    disp('excellent!')
    
    % Get pre-deployment calibration date
    prompt = 'What is the pre-deployment calibration date? [YYYY MM DD]\n';
    dnum_preC = datenum(input(prompt));
    
    % Get post-deployment calibration date
    prompt = 'What is the post-deployment calibration date? [YYYY MM DD]\n';
    dnum_posC = datenum(input(prompt));
    
    % islope = 1.0 + (b / n) [(1 / postslope) - 1.0]
    % islope = interpolated slope
    % b = number of days between pre-cruise calibration and the cast
    % n = number of days between pre- and post-cruise calibrations
    % postslope = slope from calibration sheet
    cond0 = t.conductivity;
    b = t.time - dnum_preT;
    n = dnum_posT - dnum_preT;
    
    % Get conductivity drift
    prompt = 'What is the slope correction on the post calibration?\n';
    postslope = input(prompt);
    
    % Correct conductivity drift
    islope = 1.0 + (b / n) * ((1 / postslope) - 1.0);
    t.conductivity = cond0.*islope;
    t.salinity = sw_salt(10*t.conductivity/sw_c3515,t.temperature,...
      t.pressure);
    
    % Save calibration information
    t.Cond_calibration.slope = islope;
    t.Cond_calibration.pre_cal_date = datestr(dnum_preC);
    t.Cond_calibration.pos_cal_date = datestr(dnum_posC);
    
    if params.print == 1
      % Plot salinity calibration
      fig=figure(1);
      set(gcf,'PaperUnits','centimeters')
      xSize = 30; ySize = 21;
      xLeft = (30-xSize)/2; yTop = (21-ySize)/2;
      set(gcf,'Position',[xLeft yTop xSize*40 ySize*40])
      
      subplot(2,1,1)
      plot(t.time,t.temperature,'k',t.time,temp0,'b')
      legend('Calibrated','Level 0')
      title(sprintf('%s_%s_%s_%s_%d',mooring,inst,sernum,depdep,...
        params.qcLevel),'interpreter','None','fontsize',20)
      datetick('x','mmm-yy','keeplimits')
      ylabel('Temperature (°C)')
      grid on;box on
      
      subplot(2,1,2)
      plot(t.time,t.conductivity,'k',t.time,cond0,'b')
      legend('Calibrated','Level 0')
      datetick('x','mmm-yy','keeplimits')
      ylabel('Conductivity (S/m)')
      grid on;box on
      
      set(gcf,'units','centimeters')
      set(gcf,'papersize',[32 20])
      set(gcf,'paperposition',[0,0,32,20])
      nfop=sprintf('%s_%s_%s_%s_%d',mooring,inst,sernum,depdep,...
        params.qcLevel);
      nfoppng=sprintf('%s/%s.png',params.dirpng,nfop);
      print(gcf,'-dpng','-r200',nfoppng);
      
      nfoppdf=sprintf('%s/%s.pdf',params.dirpdf,nfop);
      print(gcf,'-dpdf','-r200',nfoppdf);
      close
    end
  end
  
  % Save data to Level 1 folder
  if params.save == 1
    fname=sprintf('%s_%s_%s_%s_%d',strrep(mooring,'_',''),inst,sernum,...
      depdep,params.qcLevel);
    fpath=sprintf('%s/%s.mat',params.outpath,fname);
    save(fpath,'t')
  end
  
  % Calibrate another instrument?
  prompt='Do you want to calibrate another instrument? (y/n)\n';
  keeprunning=input(prompt,'s');
  clc
end
