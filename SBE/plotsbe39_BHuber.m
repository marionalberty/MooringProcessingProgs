function plotsbe39(s,sn, dpth)
% plot data from sbe39 serial number sn from structure s 



xl=[datenum(2011,4,4) datenum(2013,3,23)];  % plot limits
subplot 211
plot (s.ti, s.te)
datetick('x', 'mm/yy')
set(gca, 'xlim',xl)
title(['SBE39 sn ' num2str(sn) ' ' num2str(dpth) ' m']);
ylabel('Temperature [^oC]')
subplot 212
plot(s.ti, s.pr)
datetick('x', 'mm/yy')
set(gca, 'xlim',xl)
ylabel('Pressure [decibars]')

return