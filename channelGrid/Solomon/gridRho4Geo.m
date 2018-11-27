function [rho1,rho2] = gridRho4Geo(RHO1,RHO2)
% Make rho grids only using observations from CTDs deployed at the same
% depth and with data at the same time steps.

% Set max vertical distance between sensors to consider
z_buf = 50;     %[m]

%% Extract time and z grid for each RHO, get sensor overlaps, convert to
%  intermediate format
% Extract time
time1 = RHO1.time;
time2 = RHO2.time;
% Extract z
z1 = RHO1.z;
z2 = RHO2.z;
% Get all the fields that need to be gridded
fnames = fieldnames(RHO1.intermediate);
% Convert to intermediate format
for i = 1:numel(fnames)
  eval(sprintf('rho1.%s = RHO1.intermediate.%s;',fnames{i},fnames{i}))
  eval(sprintf('rho2.%s = RHO2.intermediate.%s;',fnames{i},fnames{i}))
end
% Housekeeping
clear i1 i2 i fnames RHO1 RHO2


%% Nan out observations in each record at time steps when there are no
%  observations at the same depth on the other record
% Find longer record of the two
[n_time,ilong] = max([numel(time1) numel(time2)]);
[~,ishort] = min([numel(time1) numel(time2)]);

for iL = 1:n_time
  % Step through longer record and nan out observations that don't have a
  % time, depth match
  % Find index of concurrent time step for the shorter record
  eval(sprintf('iS = find((time%i - time%i(iL)) == 0);',ishort,ilong))
  % Make sure the shorter record has observations at time(iL)
  if ~isempty(iS)
    % For the time step, find data pairs at pressure with z_buf m of
    % eachother
    % Extract profile for each record
    eval(sprintf('presS = rho%i.pres(:,iS);',ishort))
    eval(sprintf('presL = rho%i.pres(:,iL);',ilong))
    % Initialize index showing which observations to keep
    i_kS = zeros(size(presS));
    i_kL = zeros(size(presL));
    
    % Walk through pressure from the longer record
    for k = 1:numel(i_kL)
      % Find index from other observations within z_buf
      i_match = find(abs(presS - presL(k)) < z_buf);
      if ~isempty(i_match)
        i_kL(k) = 1;
        i_kS(i_match) = 1;
      end
    end
    
    % Nan out observations without match
    eval(sprintf('rho%i.pres(~i_kS,iS) = nan;',ishort))
    eval(sprintf('rho%i.sgth(~i_kS,iS) = nan;',ishort))
    eval(sprintf('rho%i.pres(~i_kL,iL) = nan;',ilong))
    eval(sprintf('rho%i.sgth(~i_kL,iL) = nan;',ilong))
  end
end

%% Make final grids
rho1 = moorInterpZ(rho1,z1,time1,'linear');
rho2 = moorInterpZ(rho2,z2,time2,'linear');