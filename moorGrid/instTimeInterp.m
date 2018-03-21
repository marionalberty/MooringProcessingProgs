function data = instTimeInterp(data,time)
% Interpolate instrument timeseries onto the mooring's common time.

% Interpolate temperature data
data.temperature = interp1(data.time,data.temperature,time);

% Interpolate pressure data
if isfield(data,'pressure')
  data.pressure = interp1(data.time,data.pressure,time);
end

% Interpolate salinity data
if isfield(data,'salinity')
  data.salinity = interp1(data.time,data.salinity,time);
end

% Interpolate velocity data
if isfield(data,'u')
  data.u = interp1(data.time,data.u',time)';
  data.v = interp1(data.time,data.v',time)';
  if isfield(data,'depth') && ~isvector(data.u)
    data.depth = interp1(data.time,data.depth',time)';
    data.pitch = interp1(data.time,data.pitch,time);
    data.roll = interp1(data.time,data.roll,time);
    data.heading = interp1(data.time,data.heading,time);
  end
end

data.time = time;