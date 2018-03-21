function [ASV,XSV] = channelInterpSS(ASV,XSV,xdist,i_data,method)
% Interpolate the XSV and ASV data between the channel sidewalls and
% mooring locations to fill the gaps using the prescribed method.
% Ignore depths (leave 0) when no velocity data is present.
% This function is specifically writen for the configuration of the
% moorings in Solomon Strait.

% Get size of velocity matrix
sz = size(ASV);

% Loop through time
for i = 1:sz(3)
  % Loop through z
  for k = 1:127
    % Extract velocity x-sections
    asv_t = ASV(k,:,i);   xsv_t = XSV(k,:,i);
    
    % Find grid points nearest to walls at this depth
    inan = diff(isnan(asv_t));
    i_start = find(inan == -1)+1;
    i_end = find(inan == 1);
    clear inan
    
    
    %% Interpolate along the first line
    
    % Use grid points near wall and between 2nd mooring and west wall
    if i_end(1) < i_data(2)
      % Topography obstruct this line
      i_start1 = i_start;
      i_end1 = [i_end(1) i_data(2)];
    else
      i_start1 = i_start(1);
      i_end1 = i_data(2);
    end
    
    for j = 1:numel(i_start1)
      % Determine what data is along the unobstructed line of depth
      i_use = i_data(i_data > i_start1(j) & i_data <= i_end1(j));
      
      % Concatenate interp points, remove duplicates
      i_all = unique([i_start1(j) i_use i_end1(j)]);
      
      if numel(i_all) > 1 && ~isempty(i_use)
        % Ensure there is sufficient data to interpolate
        
        if strcmp(method,'nearest') && numel(i_use) == 2
          % Use the data from the nearest mooring with two moorings
          asv_t(i_all(1)+1:i_all(end)) = interp1(xdist(i_use),...
            asv_t(i_use),xdist(i_all(1)+1:i_all(end)),method,'extrap');
          xsv_t(i_all(1)+1:i_all(end)) = interp1(xdist(i_use),...
            xsv_t(i_use),xdist(i_all(1)+1:i_all(end)),method,'extrap');
          
        elseif strcmp(method,'nearest') && numel(i_use) == 1
          % Use the data from the nearest mooring when there is only one
          asv_t(i_all(1)+1:i_all(end)-1) = asv_t(i_use);
          xsv_t(i_all(1)+1:i_all(end)-1) = xsv_t(i_use);
          
        elseif strcmp(method,'linear')
          % Use a linear interpolation
          asv_t(i_start1(j):i_end1(j)) = interp1(xdist(i_all),...
            asv_t(i_all),xdist(i_start1(j):i_end1(j)));
          xsv_t(i_start1(j):i_end1(j)) = interp1(xdist(i_all),...
            xsv_t(i_all),xdist(i_start1(j):i_end1(j)));
        end
      end
    end
    
    
    %% Interpolate along the second line
    
    % Update i_start and i_end to ignore points along first line
    i_start(i_start <= i_data(3)) = [];
    i_end(i_end <= i_data(3)) = [];
    
    % Append i_start if it now has fewer elements than i_end
    if numel(i_start) < numel(i_end)
      i_start = [i_data(3) i_start];
    end
    
    for j = 1:numel(i_start)
      % Determine what data is along the unobstructed line of depth
      i_use = i_data(i_data >= i_start(j) & i_data < i_end(j));
      
      % Concatenate interp points, remove duplicates
      i_all = unique([i_start(j) i_use i_end(j)]);
      
      if numel(i_all) > 1 && ~isempty(i_use)
        % Ensure there is sufficient data to interpolate
        
        if strcmp(method,'nearest') && numel(i_use) == 2
          % Use the data from the nearest mooring with two moorings
          asv_t(i_all(1)+1:i_all(end)) = interp1(xdist(i_use),...
            asv_t(i_use),xdist(i_all(1)+1:i_all(end)),method,'extrap');
          xsv_t(i_all(1)+1:i_all(end)) = interp1(xdist(i_use),...
            xsv_t(i_use),xdist(i_all(1)+1:i_all(end)),method,'extrap');
          
        elseif strcmp(method,'nearest') && numel(i_use) == 1
          % Use the data from the nearest mooring when there is only one
          asv_t(i_all(1)+1:i_all(end)-1) = asv_t(i_use);
          xsv_t(i_all(1)+1:i_all(end)-1) = xsv_t(i_use);
          
        elseif strcmp(method,'linear')
          % Use a linear interpolation
          asv_t(i_start(j):i_end(j)) = interp1(xdist(i_all),...
            asv_t(i_all),xdist(i_start(j):i_end(j)));
          xsv_t(i_start(j):i_end(j)) = interp1(xdist(i_all),...
            xsv_t(i_all),xdist(i_start(j):i_end(j)));
        end
      end
    end
    
    % Insert data back into ASX and XSV
    ASV(k,:,i) = asv_t;   XSV(k,:,i) = xsv_t;
    
  end
end