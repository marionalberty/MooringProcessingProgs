function [ASV,XSV] = channelInterpVS(ASV,XSV,xdist,i_data,method)
% Interpolate the XSV and ASV data between the channel sidewalls and
% mooring locations to fill the gaps using the prescribed method. If there
% are insufficent points for the cubic method, resort to linear
% interpolation. Ignore depths (leave 0) when no velocity data is present.
% This function is specifically writen for the configuration of the
% moorings in Vitiaz Strait.

% Get size of velocity matrix
sz = size(ASV);

% Loop through time
for i = 1:sz(3)
  % Loop through z
  for k = 1:sz(1)
    % Extract velocity x-sections
    asv_t = ASV(k,:,i);   xsv_t = XSV(k,:,i);
    
    % Find grid points nearest to walls at this depth
    inan = diff(isnan(asv_t));
    i_start = find(inan == -1)+1;
    i_end = find(inan == 1);
    clear inan
    if isempty(i_end)
      i_end = find(asv_t == 0 ,1,'last');
    end
    if isempty(i_start)
      i_start = find(asv_t == 0 ,1,'first');
    end
    
    for j = 1:numel(i_start)
      % Determine what data is along the unobstructed line of depth
      i_use = i_data(i_data > i_start(j) & i_data <= i_end(j));
      
      % Concatenate interp points, remove duplicates
      i_all = unique([i_start(j) i_use i_end(j)]);
      
      if numel(i_all) > 1 && ~isempty(i_use)
        % Ensure there is sufficient data to interpolate
        if strcmp(method,'nearest') && numel(i_use) > 1
          % Use the data from the nearest mooring with mulitple moorings
          asv_t(i_start(j)+1:i_end(j)-1) = interp1(xdist(i_use),...
            asv_t(i_use),xdist(i_start(j)+1:i_end(j)-1),method,'extrap');
          xsv_t(i_start(j)+1:i_end(j)-1) = interp1(xdist(i_use),...
            xsv_t(i_use),xdist(i_start(j)+1:i_end(j)-1),method,'extrap');
        elseif strcmp(method,'nearest') && numel(i_use) == 1
          % Use the data from the nearest mooring
          asv_t(i_start(j)+1:i_end(j)-1) = asv_t(i_use);
          xsv_t(i_start(j)+1:i_end(j)-1) = xsv_t(i_use);
        elseif strcmp(method,'linear') || ...
            (numel(i_all) < 4 && strcmp(method,'pchip'))
          % Use a linear interpolation
          asv_t(i_start(j):i_end(j)) = interp1(xdist(i_all),asv_t(i_all),...
            xdist(i_start(j):i_end(j)));
          xsv_t(i_start(j):i_end(j)) = interp1(xdist(i_all),xsv_t(i_all),...
            xdist(i_start(j):i_end(j)));
        else
          % Use a cubic interpolation
          asv_t(i_start(j):i_end(j)) = interp1(xdist(i_all),asv_t(i_all),...
            xdist(i_start(j):i_end(j)),method);
          xsv_t(i_start(j):i_end(j)) = interp1(xdist(i_all),xsv_t(i_all),...
            xdist(i_start(j):i_end(j)),method);
        end
      end
    end
    % Insert data back into ASX and XSV
    ASV(k,:,i) = asv_t;   XSV(k,:,i) = xsv_t;
  end
end