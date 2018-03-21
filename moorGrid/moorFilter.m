function DS = moorFilter(DS,time,dt)
% Filter mooring timeseries using filtfilt with a lowpass for
% frequencies lower than 1/dt. Uses a window of length dt.

% Get instrument sampling frequency
dt_inst = mean(diff(time));

% Desgin filter
w = hamming(round(dt/dt_inst));

% Get all the fields in the data structure
fnames = fieldnames(DS);

% Apply filter to all observations
for k = 1:length(fnames)
  eval(sprintf('kField = DS.%s;',fnames{k}))
  if ~isvector(kField)
    [imax,~] = size(kField);
    for i = 1:imax
      % Filter record for each instrument
      eval(sprintf('inan = ~isnan(DS.%s(i,:));',fnames{k}))
      if sum(inan)/numel(time) > 0.05
        eval(sprintf('DS.%s(i,inan) = filt_ends(w,DS.%s(i,inan));',...
          fnames{k},fnames{k}))
      end
    end
  end
end