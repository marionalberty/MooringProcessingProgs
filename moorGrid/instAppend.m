function [TT,SS,UV] = instAppend(data,TT,SS,UV)
% Append instrument timeseries onto intermediate mooring grids

%% Append intermediate temperature grid
TT.temp = [TT.temp; data.temperature(:)'];
TT.pdep = [TT.pdep; data.plannedMeterDepth];
TT.snum{end+1,1} = data.serialNo;
TT.inst{end+1,1} = data.meterType;

% Append pressure data
if isfield(data,'pressure')
  TT.pres = [TT.pres; data.pressure(:)'];
else
  TT.pres = [TT.pres; nan(size(data.temperature(:)'))];
end


%% Append intermediate salinity grid if applicable

if isfield(data,'salinity')
  SS.psal = [SS.psal; medfilt1(data.salinity(:),'omitnan')'];
  SS.pdep = [SS.pdep; data.plannedMeterDepth];
  SS.snum{end+1,1} = data.serialNo;
  SS.inst{end+1,1} = data.meterType;
  SS.pres = [SS.pres; data.pressure(:)'];
  % elseif isfield(data,'salinity_fromCTDref') && ~isfield(data,'salinity')
  %   SS.psal = [SS.psal; medfilt1(data.salinity_fromCTDref(:),'omitnan')'];
  %   SS.pdep = [SS.pdep; data.plannedMeterDepth];
  %   SS.snum{end+1,1} = data.serialNo;
  %   SS.inst{end+1,1} = data.meterType;
  %   SS.pres = [SS.pres; data.pressure(:)'];
end


%% Append intermediate velocity grid if applicable

if isfield(data,'u')
  if isvector(data.u)
    UV.U = [UV.U; data.u(:)'];
    UV.V = [UV.V; data.v(:)'];
    UV.pres = [UV.pres; data.pressure(:)'];
    UV.pdep = [UV.pdep; data.plannedMeterDepth];
    UV.snum{end+1,1} = data.serialNo;
    UV.inst{end+1,1} = data.meterType;
  else
    UV.U = [UV.U; data.u];
    UV.V = [UV.V; data.v];
    UV.pres = [UV.pres; data.depth];
    if strcmp(data.orientation,'up')
      UV.pdep = [UV.pdep; data.plannedMeterDepth-data.range];
    else
      UV.pdep = [UV.pdep; data.plannedMeterDepth+data.range];
    end
    for i=1:length(data.range)
      UV.snum{end+1,1} = data.serialNo;
      UV.inst{end+1,1} = data.meterType;
    end
  end
end