function DS = moorFillPres(DS)
% Interpolate pressure for instruments that have data but no pressure
% sensor using the planned meter depths and pressure observations of the
% other instruments on the mooring. When multiple instruments were placed
% at the same planned depth, only use one of the pressure and temperature
% records, chosing the ones with the longer timeseries.

% Buffer for removing nearby sensors
zbuf = 5;   % [m]

%% Find indecies of sensors deployed at similar depths (within 5 meters)

% Initialize indecies for which pressure and temperature sensors to use
iUse = true(size(DS.pdep));

for i = 1:length(iUse)
  % Find the vertical distance between sensors
  zdif = abs(DS.pdep(i)-DS.pdep);
  % See if any are within 5 m
  idup = find(zdif <= zbuf);
  % Check if any other planned depths are within 5 m
  if numel(idup) > 1
    % See if one of the duplicated planned depths is a SBE37
    iSBE37 = strcmp(DS.inst(idup),'SBE37');
    if sum(iSBE37) > 0
      iUse(idup) = iSBE37;
    else
      % If not, use follow algorithm
      % Count up number of nans in pressure
      pcount = sum(isnan(DS.pres(idup,:)),2);
      % Sort count of nans
      [~,ikeepP] = sort(pcount);
      % Get indecies for duplicates
      idupP = idup(ikeepP);
      % Remove signals with fewest nans
      iUse(idupP(2:end)) = false;
    end
  end
end

%% Remove unnecessary data

DS.pres = DS.pres(iUse,:);
DS.pdep = DS.pdep(iUse);
DS.snum = DS.snum(iUse);
DS.inst = DS.inst(iUse);
DS.temp = DS.temp(iUse,:);


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