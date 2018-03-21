% T_P_spectra

% Calcualte and plot spectra from sensors measuring temperature and
% pressure to determine what varibility may be due to mooring drawdown vs.
% natural variability of temperature

clear all; close all; clc

set(0,'defaultaxesfontsize',12,'defaultaxeslinewidth',0.7,...
  'defaultlinelinewidth',1,'defaultpatchlinewidth',0.7,...
  'defaultFigureColor','white')

%% Set script options

params.print=1; % 1= save figures
params.qcLevel=0; % qcLevel of data to be used

%% Set paths and directories

% Personnal paths:
% driveName='/Users/cyrilgermineaud/Documents/MATLAB/';
% addpath([driveName 'Routines_Cyril'])
driveName='/Users/marionsofiaalberty/MATLAB/Solomon_Sea/';

% Set input/output paths
params.inpath =[driveName 'Moorings/Data/Level_' num2str(params.qcLevel) '/'];

%% Calculate and plot spectra for every T-P sensor

% Get list of moorings
mlist=dir(params.inpath);

% Remove non directories
rmi=[];
for i=1:numel(mlist)
  if mlist(i).isdir ~= 1 || length(mlist(i).name) < 3
    rmi=[rmi i];
  end
end
mlist(rmi)=[];

% For each mooring, find all the T-P sensors
for i=1:numel(mlist)
  % Get mooring name
  mooring = mlist(i).name;
  
  % Set figure paths
  params.pathfig = [driveName 'Moorings/Figures/Level_' ...
    num2str(params.qcLevel) '/' mooring '/'];
  params.dirpng = [params.pathfig 'PNG'];
  params.dirpdf = [params.pathfig 'PDF'];
  
  % Get list of instruments
  ilist = dir(fullfile([params.inpath mooring '/'],'*.mat'));
  
  % Check each instrument for t & p
  for k=1:numel(ilist)
    % Load data
    load([params.inpath mooring '/' ilist(k).name])
    if exist('t','var')
      data=t;
    else
      data=s;
    end
    
    % Check for pressure
    if ~isfield(data,'pressure')
      % Go to next instrument if none
      continue
    end
    
    % Get sampling frequency [1/days]
    fs = 1./mean(diff(data.time));
    nfft = round(90*fs);
    
    % Get temp spectra
    [Tpow,freq] = fast_psd(data.temperature,nfft,fs);
    
    % Get temp spectra
    [Ppow,~] = fast_psd(data.pressure,nfft,fs);
    
    
    loglog(freq,Ppow./sum(Ppow),'b',freq,Tpow./sum(Tpow),'r')
    axis tight
    hold on
    loglog([24/12.4206 24/12.4206],[1e-7 1],'k')
    loglog([24/23.9344 24/23.9344],[1e-7 1],'k')
    loglog([24/25.8194 24/25.8194],[1e-7 1],'k')
    legend('pressure','temperature')
    xlabel('Frequency [1/day]')
    close
    
    % Clear old variables
    clear t s data
  end
  
end





