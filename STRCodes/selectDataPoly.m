function [PosInd,NegInd] = selectDataPoly(dates, temp)

% selectDataPoly: from Ron's code selectSBE39

% USAGE: within cleanSTR code
% Oliver Vetter: November 2012

% Begin selecing points, using a n-point polygon;
% Initially, the list of points is empty.
x = [];
y = [];
n = 0;
but = 1;
while but == 1
    [xi,yi,but] = ginput(1);
    plot(xi,yi,'+r');
    n = n+1;
    x(n) = xi;
    y(n) = yi;
end
ind = inpolygon(dates, temp, x,y);

% Positive Indicies - data NOT selected for removal
PosInd = find(ind == 0);  
% Negative Indicies - data seclected for removal
NegInd = find(ind == 1);  

% plots the removable selected data
plot(dates(NegInd),temp(NegInd),'-+r');