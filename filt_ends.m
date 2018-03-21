function fx = filt_ends(wi,x)
% function fx = filt_ends(wi,x)
%
% filters using filtfilt but does a better job by
% buffering data at two ends:
% Buffer method:  
%   bufmode = 0 - mirrors data 2*length(wi) at each end [DEFAULT]
%   bufmode = -1 - turns off buffering
% NOTE: can't tolerate NaNs
%
% Susan Wijffels
% CMR
% Sept 19, 2001
lf = length(wi);
nx= length(x);
frontend = x(2*lf:-1:1);
backend = x([nx-2*lf]:nx);
xx = [frontend(:);x(:);backend(:)];
keep = find(~isnan([NaN*frontend(:);x(:);NaN*backend(:)]));
xxf = filtfilt(wi,sum(wi),xx);
fx = xxf(keep);
return
