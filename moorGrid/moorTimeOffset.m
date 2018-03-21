function DS=moorTimeOffset(DS,useField,time)
% Correct for the temporal offset between instruments along the mooring by
% working on adjacent pairs, assuming the first sensor is correct and that
% adjacent sensors move roughly in phase with tide induced motions.
% Only works for signals with the same sampling frequency.
% DS = data structure (either for temp or vel)
% useField = name of field to use for correlation

% Set sampling frequency
FS=1/3600;    % [Hz]

% Initialize time offset
offset = zeros(size(DS.pdep));

% Get offset between pairs
for i = 2:numel(offset)
  % Get index from when both instruments have started recording
  eval(sprintf('i1 = find(~isnan(DS.%s(i-1,:)),1);',useField))
  eval(sprintf('i2 = find(~isnan(DS.%s(i,:)),1);',useField))
  istart = max([i1 i2]);
  % Get index for 90 days later
  iend = istart+(90*24);
  
  % Get signals
  eval(sprintf('signal1 = DS.%s(i-1,istart:iend);',useField))
  eval(sprintf('signal2 = DS.%s(i,istart:iend);',useField))
  
  % Estimate the coherence over the desired time period
  % nfft = 256; (greater than 1 week and 2^8 for fft speed)
  p=fast_cohere(signal1,signal2,256,FS);
  % Find the low frequencies for which the signals are coherent
  ginds=find(p.coh>.5 & p.f<0.2);
  % Unwrap phase for total time offset at each frequency
  temp=unwrap(p.pha);
  % Fit a line to only the coherent offsets
  b=regress(temp(ginds)',p.f(ginds)');
  
  % Get offset back into days
  offset(i)=b(1)/2/pi/24/3600;
end

% Get cumlative offset compared to grid time
total_offset=cumsum(offset);

% Get all the fields that need to be corrected
fname = fieldnames(DS);

% Correct each field of the data structure
for k = 1:numel(fname)
  % Ensure not a vector
  eval(sprintf('field2fix = DS.%s;',fname{k}))
  
  if ~isvector(field2fix)
    % Correct each timeseries within the field for it's offset
    for i = 2:numel(offset)
      eval(sprintf(...
        'DS.%s(i,:) = interp1(time+total_offset(i),DS.%s(i,:),time,''linear'',''extrap'');',...
        fname{k},fname{k}))
    end
  end
end