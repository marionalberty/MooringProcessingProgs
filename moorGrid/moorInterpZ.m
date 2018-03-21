function DS_out = moorInterpZ(DS,z,time,method)
% Interpolate the final version of the intermediate level data onto
% regularly spaced z  and time grid.

% Make DS_out have a structure that contains the intermediate level of
% the data for future reference.
DS_out.intermediate = DS;
% Add time and z
DS_out.z = z;
DS_out.time = time;


%% Initialize final data grids

% Get all the fields that need to be gridded
fname = fieldnames(DS);
% Remove other fields from list
fname(strcmp(fname,'pres')) = [];
fname(strcmp(fname,'snum')) = [];
fname(strcmp(fname,'inst')) = [];
fname(strcmp(fname,'pdep')) = [];

for i = 1:numel(fname)
  eval(sprintf('DS_out.%s = nan(length(z),length(time));',fname{i}))
end


%% Interp (but don't extrap) each data grid onto regular z grid

% More than one inst to work with
for i = 1:numel(fname)
  for k = 1:length(time)
    % Extract pressure vector
    p = DS.pres(:,k);
    % Get indecies for data to use in the interpolation that has a
    % corresponding pressure record
    eval(sprintf('i_data = find(~isnan(p) & ~isnan(DS.%s(:,k)));',...
      fname{i}))
    if numel(i_data) > 1
      % Interp data
      if strcmp(method,'pchip')
        eval(sprintf(...
          'DS_out.%s(:,k) = interp1(p(i_data),DS.%s(i_data,k),z,''pchip'',nan);',...
          fname{i},fname{i}))
      elseif strcmp(method,'linear')
        eval(sprintf('DS_out.%s(:,k) = interp1(p(i_data),DS.%s(i_data,k),z);',...
          fname{i},fname{i}))
      end
    elseif numel(i_data) == 1
      % Get closest z bin
      [~,i_zmin] = min(abs(z-p(i_data))); %#ok<ASGLU>
      % Assign the variable's value to the nearest bin
      eval(sprintf('DS_out.%s(i_zmin,k) = DS.%s(i_data,k);',fname{i},...
        fname{i}))
    end
  end
end