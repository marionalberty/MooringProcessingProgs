function [powe,fre]=fast_psd(x,nfft_in,fs)
% function [p,f]=fast_psd(x,nfft,Fs) computes the properly
% normalized psd for a series x.  Just like matlab's.  nfft is the
% desired fourier transform length, and Fs is the sampling
% frequency.  This routine detrends the data and applies a hanning
% window before transforming.  A minium of 50% overlap is used;
% segments are fit so that all of the data is used.

% ensure nfft_in is greater than segment length
max_ind=length(x);
nfft=min([nfft_in max_ind]);
% number of segments with 50% overlap
repeats=fix(2*max_ind/nfft);
if max_ind==nfft
  repeats=1;
  if nfft~=nfft_in
    % warning(['Only ' num2str(nfft) ' points are being used for this spectra']);
  end
end
% create hannging window
wind=hanning(nfft);
if size(x,1)==1
  wind=wind';
end


% I believe the following is the correct normalization for the PSD.
W1 = 2/norm(wind)^2 ;
% jonathan's psd routine which is about twice the speed of the canned
% one....

% preform fft on first segment
total=fft(detrend(x(1:nfft)).*wind);
% get power from first fft segment
powe=total(2:floor(nfft/2+1)).*conj(total(2:floor(nfft/2+1)));
% repeat for remaining segments
if (repeats-1)
  step=fix((max_ind-nfft)/(repeats-1));
  for i=step:step:(max_ind-nfft)
    % compute fft for next segment
    total=fft(detrend(x(i:(i+nfft-1))).*wind);
    % add power from above computed segment to the power of all
    % previous segments
    powe=powe+total(2:nfft/2+1).*conj(total(2:nfft/2+1));
  end
end
% multiple by 2 to account for only using half of the fft, divide by
% repeat to get mean power for all segments and normalize by sampling
% frequency
powe = W1*powe/repeats/fs;
% get appropriate frequency range
fre = linspace(fs/nfft,fs/2,nfft/2)';
if size(x,1)==1
  fre=fre';
end