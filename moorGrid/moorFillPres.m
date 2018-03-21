function DS = moorFillPres(DS)
% Interpolate pressure for instruments that have data but no pressure
% sensor using the planned meter depths and pressure observations of the
% other instruments on the mooring. When multiple instruments were placed
% at the same planned depth, only use one of the pressure and temperature
% records, chosing the ones with the longer timeseries.

% Buffer for removign nearby sensors
zbuf = 5;   % [m]

%% Find indecies of sensors deployed at similar depths (within 5 meters)

% Initialize indecies for which pressure and temperatyre sensors to use
iUseP = true(size(DS.pdep));
iUseT = true(size(DS.pdep));

for i = 1:length(iUseP)
  zdif = abs(DS.pdep(i)-DS.pdep);
  idup = find(zdif <= zbuf);
  if ~isscalar(idup)
    tcount = sum(isnan(DS.temp(idup,:)),2);
    pcount = sum(isnan(DS.pres(idup,:)),2);
    [~,ikeepT] = sort(tcount);
    [~,ikeepP] = sort(pcount);
    idupT = idup(ikeepT);
    idupP = idup(ikeepP);
    iUseT(idupT(2:end)) = false;
    iUseP(idupP(2:end)) = false;
  end
end

%% Remove unnecessary data

DS.pres = DS.pres(iUseP,:);
DS.pdep = DS.pdep(iUseP);
DS.snum = DS.snum(iUseP);
DS.inst = DS.inst(iUseP);
DS.temp = DS.temp(iUseT,:);


%% Interpolate for pressure at each timestep for instruments with temp data
%  but not pressure.

for i = 1:length(DS.pres(1,:))
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