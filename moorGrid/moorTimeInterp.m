function DS = moorTimeInterp(DS,time_initial,time_final)
% Interpolate mooring timeseries onto a new (typically less frequent) time
% vector.

% Get all the fields in the data structure
fnames = fieldnames(DS);

% Apply interpolation to all observations
for k = 1:length(fnames)
  eval(sprintf('kField = DS.%s;',fnames{k}))
  if ~iscolumn(kField)
    eval(sprintf('DS.%s = interp1(time_initial,DS.%s'',time_final)'';',...
      fnames{k},fnames{k}))
    eval(sprintf('kField = DS.%s;',fnames{k}))
    if iscolumn(kField)
       eval(sprintf('DS.%s = DS.%s'';',fnames{k},fnames{k}))
    end
  end
end