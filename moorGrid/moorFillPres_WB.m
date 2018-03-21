function DS = moorFillPres_WB(DS,params,time)
% SAME AS moorFillPres.m BUT HANDLES A MOORING THAT EXPERIENCED A LINE
% BREAK.

% Interpolate pressure for instruments that have data but no pressure
% sensor using the planned meter depths and pressure observations of the
% other instruments on the mooring. When multiple instruments were placed
% at the same planned depth, only use one of the pressure and temperature
% records, chosing the ones with the longer timeseries. Additionally,
% interpolate per usual the timeseries prior to the break but only
% interpolate pressure for instruments below the fall depth after the break
% occured.

% Buffer for removign nearby sensors
zbuf = 5;   % [m]
% Find the indecies just before and after the break
iprirB = find(time < params.breaktime,1,'last');
ipostB = find(time > params.breaktime,1);

%% Find indecies of sensors deployed at similar depths (within 5 meters)

% Initialize indecies for which pressure and temperatyre sensors to use
iUseP = true(size(DS.pdep));
iUseT = true(size(DS.pdep));

for i = 1:length(iUseP)
  % Get difference in planned meter depth
  zdif = abs(DS.pdep(i)-DS.pdep);
  % Find differences less than the buffer
  idup = find(zdif <= zbuf);
  if ~isscalar(idup)
    % If idup is a vector, find which sensors has more observations 
    tcount = sum(isnan(DS.temp(idup,:)),2);
    pcount = sum(isnan(DS.pres(idup,:)),2);
    [~,ikeepT] = sort(tcount);
    [~,ikeepP] = sort(pcount);
    idupT = idup(ikeepT);
    idupP = idup(ikeepP);
    % Ignore data with worse coverage
    iUseT(idupT(2:end)) = false;
    iUseP(idupP(2:end)) = false;
  end
end

%% Remove unnecessary data

DS.pres = DS.pres(iUseP,:);
DS.pdep = DS.pdep(iUseT);
DS.snum = DS.snum(iUseT);
DS.inst = DS.inst(iUseT);
DS.temp = DS.temp(iUseT,:);

% Get indecies of instruments that recorded temp but not pres
i_nopr = find(sum(isnan(DS.pres),2)./length(time) == 1);


%% Interpolate for pressure at each timestep for instruments with temp data
%  but not pressure for times prior to the mooring break

for i = 1:iprirB
  % Extract column vector
  p = DS.pres(:,i);
  t = DS.temp(:,i);
  % Get indecies for instruments with temp data
  i_wtemp = find(~isnan(t));
  % Get indecies for instruments with just pressure data
  i_wpres = find(~isnan(p));
  % Interpolate using planned meter depths
  DS.pres(i_wtemp,i) = interp1(DS.pdep(i_wpres),p(i_wpres),...
    DS.pdep(i_wtemp),'linear','extrap');
end


%% Interpolate for pressure at each timestep for instruments with temp data
%  but not pressure for times after the mooring break

% Get index for instruments that fell
ifall = find(DS.pdep < params.fallDepth);
% get index of fallen instrutments w/o pressure
i_nanT = find(ismember(ifall,i_nopr));
% Nan out temperature observations from fallen instruments without pressure
% sensors.
DS.temp(i_nanT,ipostB:end) = nan;

for i = ipostB:length(time)
  % Extract column vector
  p = DS.pres(:,i);
  t = DS.temp(:,i);
  % Get indecies for instruments with temp data
  i_wtemp = find(~isnan(t));
  % Get indecies for instruments with just pressure data
  i_wpres = find(~isnan(p));
  % Interpolate using planned meter depths
  DS.pres(i_wtemp,i) = interp1(DS.pdep(i_wpres),p(i_wpres),...
    DS.pdep(i_wtemp),'linear','extrap');
end