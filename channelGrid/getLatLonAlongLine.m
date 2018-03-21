function [x_out,y_out] = getLatLonAlongLine(x,y,mb)

% x  = longitude vector
% y  = latitude vector
% mb = [slope intercept] of desired line

% Line equations: y = m*x+b
line_y = @(x,mb) mb(1).*x + mb(2);
line_x = @(y,mb) (y-mb(2))./mb(1);

% Calc lines for available x and y grid
y_trueX = line_y(x,mb);
x_trueY = line_x(y,mb);

% Find the nearest available grid point for each version of the line
y_coordX = nan(size(x));
for i = 1:numel(x)
  [~,i_y] = min(abs(y_trueX(i) - y));
  y_coordX(i) = y(i_y);
end

x_coordY = nan(size(y));
for i = 1:numel(y)
  [~,i_x] = min(abs(x_trueY(i) - x));
  x_coordY(i) = x(i_x);
end

% Concatinate the list of x,y coordinates
xy_out = [x' y_coordX'; x_coordY' y'];

% Remove duplicate coordinate pairs
xy_out = unique(xy_out,'rows');
x_out = xy_out(:,1)';
y_out = xy_out(:,2)';