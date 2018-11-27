% extrap_mooringTimeseries_GCUR.m

% Extrapolate mooring velocity observations to the surface and the sea
% floor for the Solomon Strait moorings utilizing GlobCurrent total surface
% velocity (geo + ekman) to constrain surface extrapolation.

% Extrapolate mooring velocity observations to the surface and the sea
% floor.

clear all; close all; clc
warning off

set(0,'defaultaxesfontsize',16,'defaultaxeslinewidth',0.7,...
  'defaultlinelinewidth',1,'defaultpatchlinewidth',0.7,...
  'defaultFigureColor','white')

%% Set paths and directories
% Personnal paths:
% driveName='/Users/cyrilgermineaud/Documents/MATLAB/';
% addpath([driveName 'Routines_Cyril'])
driveName = '/Users/marionsofiaalberty/MATLAB/';
dataPath = 'Solomon_Sea/Moorings/Data/Gridded/';
% switch for plottting
i_plot = 0;
% Surface extrap
out_surf = 'gcur4surf';


% Choose cases
for i_moor = 1:3
  for i_bottom = 1:2
    %% Set up each mooring case and load data
    switch i_moor
      case 1
        % Solomon M1
        params.moor_name = 'Solomon_M1';
        params.channel = 'SolomonStrait';
        load([driveName dataPath params.channel '/' params.moor_name ...
          '/Solomon_M1.mat'])
      case 2
        % Solomon M2b
        params.moor_name = 'Solomon_M2b';
        params.channel = 'SolomonStrait';
        load([driveName dataPath params.channel '/' params.moor_name ...
          '/Solomon_M2b.mat'])
        load([driveName dataPath params.channel '/Solomon_M2a' ...
          '/Solomon_M2a.mat'],'SS','TT','RHO','T4S')
      case 3
        % Solomon M3
        params.moor_name = 'Solomon_M3';
        params.channel = 'SolomonStrait';
        load([driveName dataPath params.channel '/' params.moor_name ...
          '/Solomon_M3.mat'])
    end
    
    %% Add and adjust params
    % Move the intermediate UV data (pre-extrapolation data)
    params.intermediateUV = UV.intermediate;
    UV = rmfield(UV,'intermediate');
    
    
    %% Load velocities
    load([driveName dataPath params.channel '/' params.moor_name ...
      '/UV_GlobCurrents.mat'],'UV_SC')
    
    
    %% Extrap velocity to surface
    % Interpolate from the shallowest velocity to surface velocity
    for i = 1:length(UV.time)
      % Find the end of the nans near the surface
      i_top = find(diff(isnan(UV.U(:,i))) == -1);
      if ~isempty(i_top)
        % Linearly interpolate between the surface and the shallowest
        % velocity
        UV.U(1:i_top,i) = interp1([0; UV.z(i_top+1)],...
          [UV_SC.U(i); UV.U(i_top+1,i)],UV.z(1:i_top));
        UV.V(1:i_top,i) = interp1([0; UV.z(i_top+1)],...
          [UV_SC.V(i); UV.V(i_top+1,i)],UV.z(1:i_top));
      end
    end
    
    %% Save surface Ekman and Geostrophic currents
    UV.U_geo = UV_SC.U_geo;
    UV.V_geo = UV_SC.V_geo;
    UV.U_ekm = UV_SC.U_ekm;
    UV.V_ekm = UV_SC.V_ekm;
    
    
    %% Extrap velocity to the bottom
    % Since the last grid point is either the bottom or sill depth,
    % set the velocity there to zero
    UV.U(end,:) = 0;
    UV.V(end,:) = 0;
    
    switch i_bottom
      case 1
        % Constant velocity with zero velocity in last bin
        out_bot = 'slab2bot';
        for i = 1:length(UV.time)
          % Find the start of the nans near the bottom
          i_bot = find(diff(isnan(UV.U(:,i))) == 1)+1;
          % Apply slab
          UV.U(i_bot:end-1,i) = UV.U(i_bot-1,i);
          UV.V(i_bot:end-1,i) = UV.V(i_bot-1,i);
        end
      case 2
        % Constant shear to zero
        out_bot = 'linear2bot';
        for i = 1:length(UV.time)
          % Find the start of the nans near the bottom
          i_bot = find(diff(isnan(UV.U(:,i))) == 1)+1;
          % Linearly interpolate to zero
          UV.U(i_bot-1:end,i) = interp1([UV.z(i_bot-1); UV.z(end)],...
            [UV.U(i_bot-1,i); UV.U(end,i)], UV.z(i_bot-1:end));
          UV.V(i_bot-1:end,i) = interp1([UV.z(i_bot-1); UV.z(end)],...
            [UV.V(i_bot-1,i); UV.V(end,i)], UV.z(i_bot-1:end));
        end
    end
    
    
    %% Add magnitude, direction, ASV, and XSV to structure
    % Magnitude
    UV.magnitude = abs(UV.U + 1i*UV.V);
    % Direction (cardinal)
    UV.direction = wrapTo360(rad2deg(atan2(UV.U,UV.V)));
    
    
    %% Save data
    % finalize params
    params.bottom_method = out_bot;
    params.surface_method = out_surf;
    params.print = i_plot;
    
    fname = [driveName dataPath params.channel '/' params.moor_name ...
      '/' params.moor_name '_' out_bot '_' out_surf '.mat'];
    save(fname,'UV','params','SS','TT','RHO','T4S')
    
    
    %% Plot data
    if params.print
      moorPlotFinal(UV,SS,TT,params)
    end
    close all
  end
end