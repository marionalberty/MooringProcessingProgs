function [ASV,XSV] = channelInterpSTG(ASV,XSV,xdist,i_data,method)
% Interpolate the XSV and ASV data between the channel sidewalls and
% mooring locations to fill the gaps using the prescribed method. If there
% are insufficent points for the cubic method, resort to linear
% interpolation. Ignore depths (leave 0) when no velocity data is present.
% This function is specifically writen for the configuration of the 
% moorings in St George's Channel.

% Get size of velocity matrix
sz = size(ASV);

% Loop through time
for i = 1:sz(3)
  % Loop through z
  for k = 1:sz(1)
    % Extract velocity x-sections
    asv_t = ASV(k,:,i);   xsv_t = XSV(k,:,i);
    
    % Find grid point nearest to wall at that depth
    i_start = find(asv_t == 0,1,'first');
    i_end = find(asv_t == 0 ,1,'last');
    
    % Determine what data is inbetween i_start and i_end
    i_use = i_data(i_data > i_start & i_data < i_end);
    
    % Concatenate interp points, remove duplicates    
    i_all = unique([i_start i_use i_end]);
    
    % Move on to next depth if i_all has less than three points
    if numel(i_all) < 3
      continue
    end
    
    if strcmp(method,'nearest') && numel(i_all) > 3
      % Use the data from the nearest mooring with two moorings
      asv_t(i_start+1:i_end-1) = interp1(xdist(i_data),...
        asv_t(i_data),xdist(i_start+1:i_end-1),method,'extrap');
      xsv_t(i_start+1:i_end-1) = interp1(xdist(i_data),...
        xsv_t(i_data),xdist(i_start+1:i_end-1),method,'extrap');
      
    elseif strcmp(method,'nearest') && numel(i_all) == 3
      % Use the data from the nearest mooring when there is only one
      asv_t(i_start+1:i_end-1) = asv_t(i_use);
      xsv_t(i_start+1:i_end-1) = xsv_t(i_use);
      
    elseif strcmp(method,'linear') || ...
        (numel(i_all) < 4 && strcmp(method,'pchip'))
      % Use a linear interpolation
      asv_t(i_start:i_end) = interp1(xdist(i_all),asv_t(i_all),...
        xdist(i_start:i_end));
      xsv_t(i_start:i_end) = interp1(xdist(i_all),xsv_t(i_all),...
        xdist(i_start:i_end));
      
    else
      % Use a cubic interpolation
      asv_t(i_start:i_end) = interp1(xdist(i_all),asv_t(i_all),...
        xdist(i_start:i_end),method);
      xsv_t(i_start:i_end) = interp1(xdist(i_all),xsv_t(i_all),...
        xdist(i_start:i_end),method);
      
    end
    % Insert data back into ASX and XSV
    ASV(k,:,i) = asv_t;   XSV(k,:,i) = xsv_t;
    
  end
end