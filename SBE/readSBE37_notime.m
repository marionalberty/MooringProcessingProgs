function s = readSBE37_notime(file, lat, lon, press, getCond)

% readSBE37 - read a Microcat (SBE37) data logger output 'asc' file.
%
% s = readSBE37(file) reads the SBE37 output 'asc' file and creates the
% output structure 's' with the following fields of length equal to the
% number of samples:
%
%   time        - the sample time (in Matlab date format).
%   temperature - the sample temperature (deg. C ITS-90).
%   salinity    - the salinity if the Microcat contains a conductivity
%                 sensor (PSU).
%   pressure    - the sample pressure if the Microcat contains a pressure
%                 sensor (db).
%   units       - cell array containing the parameter name and units for
%                 each parameter.
% s = readSBE37(file, lat, lon) in addition inserts scalar fields for the
% location of the instrument, i.e. fields:
%
% latitude      - the instrument latitude (lat - decimal deg.).
% longitude     - the instrument longitude (lat - decimal deg.).
% s = readSBE37(file, lat, lon, press) is to be used for instruments
% without a pressure gauge and adds the scalar pressure field equal to
% press (db).  This pressure is also used to determine salinity from
% conductivity (default zero).
% s = readSBE37(file, lat, lon, press, getCond) in addition returns
% conductivity in the structure if 'getCond' is true.

% Check input parameters

error(nargchk(1, 5, nargin));

if ~ischar(file)
    error('Parameter ''file'' must be a string');
end

if isempty(findstr(file, '.asc'))
    error('File is not of ''asc'' type');
end

if nargin > 1
    if ~isnumeric(lat)
        error('Invalid latitude');
    end
else
    lat = [];
end

if nargin > 2
    if ~isnumeric(lon)
        error('Invalid longitude');
    end
else
    lon = [];
end

if nargin > 3
    if ~isnumeric(press)
        error('Invalid pressure');
    end
else
    press = 0.0;
end

if nargin < 5
    getCond = 0;
end

% Open the file and read header

fid = fopen(file, 'r');
if fid < 0
    error(['Unable to open file ', file]);
end

nSensors = 0;
tempIndex = -1;
condIndex = -1;
pressIndex = -1;
form = '%f';
line = fgetl(fid);
while isempty(line) | length(line) < 5 | line(1) == '*' | line(1) == 's'
    if ~[sum(line == ' ') == length(line) ];   % do not process if all blanks  SEW Dec 2006
        if ~isempty(findstr(line, '* temperature')) & tempIndex < 0
            nSensors = nSensors + 1;
            tempIndex = nSensors;
        else
            if ~isempty(findstr(line, '* conductivity')) & condIndex < 0
                nSensors = nSensors + 1;
                condIndex = nSensors;
                form = [form, ', %f'];
            else
                if ~isempty(findstr(line, '* pressure')) & pressIndex < 0
                    nSensors = nSensors + 1;
                    pressIndex = nSensors;
                    form = [form, ', %f'];
                end
            end
         end
         % pull out start date and sample interval
         if findstr(line,'start time'), 
            timestart = datenum(line(14:end));
         end
         if findstr(line,'sample interval'),
            deltatime = sscanf(line(18:22),'%f')/[3600*24] % converted from seconds to days
         end
    end
    line = fgetl(fid);
 end
 datestr(timestart)
 deltatime*24*3600

if tempIndex < 0
    error('Invalid file format');
end

% Read and decode the data

t = [];
c = [];
p = [];
time = [];
if pressIndex < 0
    p = press;
end

n = 0;
while ischar(line)
    n = n + 1;
    if mod(n - 1, 1000) == 0
        disp(line)
        disp(['Reading scan no. ', num2str(n)]);
    end
    
    [a, count, msg, k] = sscanf(line, form, nSensors);
    if ~isempty(msg)
        error(['Invalid input line: ', line]);
    end

    % Decode the date/time and store data

% time(n) = datenum(line(k+2:end));
    time(n) = timestart + (n-1)*deltatime;
    if tempIndex > 0
        t(n) = a(tempIndex);
    end
    
    if condIndex > 0
        c(n) = a(condIndex);
    end
    
    if pressIndex > 0
        p(n) = a(pressIndex);
    end
    
    line = fgetl(fid);
end

fclose(fid);

% Create the output structure

if condIndex > 0
    salt = sw_salt(c / 4.2914, t, p);
end

s = struct('time', time, 'temperature', t);
units = {'time', 'Matlab date'; 'temperature', 'ITS-90, deg C'};
nUnits = 2;
if condIndex > 0
    s = setfield(s, 'salinity', salt);
    nUnits = nUnits + 1;
    units(nUnits, :) = {'salinity', 'PSU'};
    if getCond
        s = setfield(s, 'conductivity', c);
        nUnits = nUnits + 1;
        units(nUnits, :) = {'conductivity', 'S/m'};
    end
end

if pressIndex > 0
    s = setfield(s, 'pressure', p);
    nUnits = nUnits + 1;
    units(nUnits, :) = {'pressure', 'db'};
else
    if nargin > 3
        s = setfield(s, 'pressure', press);
        nUnits = nUnits + 1;
        units(nUnits, :) = {'pressure', 'db'};
    end
end

s=setfield(s,'units', units);
if ~isempty(lat) & ~isempty(lon)
    s = setfield(s, 'latitude', lat);
    s = setfield(s, 'longitude', lon);
end
