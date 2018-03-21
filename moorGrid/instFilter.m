function data = instFilter(data,dt)
% Filter instrument timeseries using filtfilt with a lowpass for
% frequencies lower than 1/dt. Uses a window of length dt.

% Get instrument sampling frequency
dt_inst = mean(diff(data.time));

% Desgin filter
w = hamming(round(dt/dt_inst));

% Apply filter to temperature
inan = ~isnan(data.temperature);
data.temperature(inan) = filt_ends(w,data.temperature(inan));

% Apply filter if pressure exists
if isfield(data,'pressure')
  inan = ~isnan(data.pressure);
  data.pressure(inan) = filt_ends(w,data.pressure(inan));
end

% Apply filter if salinity exists
if isfield(data,'salinity')
  inan = ~isnan(data.salinity);
  data.salinity(inan) = filt_ends(w,data.salinity(inan));
end

% Apply filter if u and v exist
if isfield(data,'u')
  inan = ~isnan(data.u);
  data.u(inan) = filt_ends(w,data.u(inan));
  data.v(inan) = filt_ends(w,data.v(inan));
end