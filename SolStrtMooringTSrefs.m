% SolomonMooringTSrefs.m

% Use CTD profiles of T and S to come up with linear relationships of T and
% S over the range of observed T for deep moorings in attempt to estimate
% density and thus geostrophic velocity.
clear all; clc; close all
warning off

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',1,...
  'defaultlinelinewidth',1,'defaultpatchlinewidth',1,...
  'defaultFigureColor','white')


%% Set Parameters
% Personnal paths:
driveName = '/Users/marionsofiaalberty/MATLAB/Solomon_Sea/';
moorPath = [driveName 'Moorings/Data/Level_1/'];
pandoPath = [driveName 'Pandora2012/CTD/procMSA/'];
mspicePath = [driveName 'moorSPICE2014/CTD/procMSA/'];


%% Set mooring case
for i_moor = 1:3
  switch i_moor
    case 1
      % Solomon M1
      mooring = 'Solomon_M1';
      % List of CTD flist
      plist = {'pn05201'};
      mlist = {'TN307_C07_01'};
      % List of mooring sensors to get relationship for
      tsens_woS = {'SolomonM1_SBE39_5631_1400_1',...
        'SolomonM1_SBE37_7969_1700_1'};
      
    case 2
      % Solomon M2a
      mooring = 'Solomon_M2a';
      % List of CTD flist
      plist = {'pn04201','pn04202','pn04204','pn04205','pn04207',...
        'pn04208','pn04209'};
      mlist = {'TN307_C04_01','TN307_C44_01'};
      % List of mooring sensors to get relationship for
      tsens_woS = {'SolomonM2a_SBE39_1105_1400_1',...
        'SolomonM2a_SBE39_1106_2000_1','SolomonM2a_SBE39_1107_2500_1'};
      
    case 3
      % Solomon M3
      mooring = 'Solomon_M3';
      % List of CTD flist
      plist = {'pn04801'};
      mlist = {'TN307_C01_01'};
      % List of mooring sensors to get relationship for
      tsens_woS = {'SolomonM3_SBE39_1104_2000_1',...
        'SolomonM3_SBE37_9119_1700_1'};
  end
  
  %% Load CTD T-S and make reference dataset
  % Initialize T and S
  TT = [];
  SS = [];
  
  % Get Pandora data
  for i_file = 1:numel(plist)
    load([pandoPath plist{i_file} '.mat'],'datad_1m')
    TT = [TT; datad_1m.t1(datad_1m.t1 < 4)];
    SS = [SS; datad_1m.s1(datad_1m.t1 < 4)];
  end
  % Get moorSPICE data
  for i_file = 1:numel(mlist)
    load([mspicePath mlist{i_file} '.mat'],'datad_1m')
    TT = [TT; datad_1m.t1(datad_1m.t1 < 4)];
    SS = [SS; datad_1m.s1(datad_1m.t1 < 4)];
  end
  % Housekeeping
  clear datad_1m
  
  
  %% Load mooring data and get T-S relationship
  for i_sens = 1:numel(tsens_woS)
    load([moorPath mooring '/' tsens_woS{i_sens} '.mat'],'t')
    % Get t_range to calculate relationship over
    t_range = [max(t.temperature) min(t.temperature)];
    % Get indecies of reference data that is within range
    i_4calc = find(TT <= t_range(1) & TT >= t_range(2));
    % Calculate the correlation coeficients and pause
    cc_TS = corrcoef([TT(i_4calc) SS(i_4calc)]);
    % Get coefficients for relationship
    p_TS = polyfit(TT(i_4calc),SS(i_4calc),1);
    % Calculate mooring salinity from CTD refs
    t.salinity_fromCTDref = polyval(p_TS,t.temperature);
    t.corrcoef_CTD = cc_TS(2);
    % Save level 1 data
    save([moorPath mooring '/' tsens_woS{i_sens} '.mat'],'t')
  end
end
