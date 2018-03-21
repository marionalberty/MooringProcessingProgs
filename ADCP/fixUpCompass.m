function adcp_UP = fixUpCompass(adcp_UP,adcp_DOWN,i_moor)
% i_moor is a switch since the pitch sensor for StGW had issues
% i_moor = 1 : StGE
% i_moor = 2 : StGW

%% Step through time and rotate all velocity into beam coordinates

% Initialize beam velocity matricies
beam1_vel = nan(size(adcp_UP.east_vel));
beam2_vel = nan(size(adcp_UP.east_vel));
beam3_vel = nan(size(adcp_UP.east_vel));

% loop through time
for i = 1:length(adcp_UP.mtime)
  % Calc inverse of transformation matrix
  % enu = T*beam
  % beam = inv(T)*enu
  T = transformADCPcoord(adcp_UP.heading(i),adcp_UP.pitch(i),...
    adcp_UP.roll(i));
  % loop through each bin at each time step
  for k = 1:adcp_UP.config.n_cells
    beam = T\[adcp_UP.east_vel(k,i); adcp_UP.north_vel(k,i); ...
      adcp_UP.vert_vel(k,i)];
    beam1_vel(k,i) = beam(1);
    beam2_vel(k,i) = beam(2);
    beam3_vel(k,i) = beam(3);
  end
end

%% Interpolate Down ADCP heading, pitch, roll onto Up time

heading_DOWN = interp1(adcp_DOWN.mtime,adcp_DOWN.heading,adcp_UP.mtime);
pitch_DOWN = interp1(adcp_DOWN.mtime,adcp_DOWN.pitch,adcp_UP.mtime);
roll_DOWN = interp1(adcp_DOWN.mtime,adcp_DOWN.roll,adcp_UP.mtime);

%% Calculate heading offset from pitch and roll

switch i_moor
  case 1
    % St George's East
    % Use all since no clear issues with pitch and roll
    hoff = fminsearch('checktilt',0,[],[adcp_UP.roll; adcp_UP.pitch; ...
      roll_DOWN; pitch_DOWN]);
  case 2
    % St George's West
    % Only use a subset since pitch maxes out for UP in the first part of
    % the deployment.
    hoff = fminsearch('checktilt',0,[],[adcp_UP.roll(9000:13500); ...
      adcp_UP.pitch(9000:13500); roll_DOWN(9000:13500); ...
      pitch_DOWN(9000:13500)]);
end

%% Step through time and rotate all velocity back to earth with a compass
%  correction added to the heading

% Initialize beam velocity matricies
east_vel_C = nan(size(beam1_vel));
north_vel_C = nan(size(beam1_vel));
vert_vel_C = nan(size(beam1_vel));

% loop through time
for i = 1:length(adcp_UP.mtime)
  % Calc inverse of transformation matrix
  % enu = T*beam
  % beam = inv(T)*enu
  T = transformADCPcoord(heading_DOWN(i)+hoff,adcp_UP.pitch(i),...
    adcp_UP.roll(i));
  % loop through each bin at each time step
  for k = 1:adcp_UP.config.n_cells
    enu = T * [beam1_vel(k,i); beam2_vel(k,i); beam3_vel(k,i)];
    east_vel_C(k,i) = enu(1);
    north_vel_C(k,i) = enu(2);
    vert_vel_C(k,i) = enu(3);
  end
end

%% Save rotated velocities to a new file
% Add new velocity and corrected heading
adcp_UP.east_vel = east_vel_C;
adcp_UP.north_vel = north_vel_C;
adcp_UP.vert_vel = vert_vel_C;
adcp_UP.heading = heading_DOWN+hoff;