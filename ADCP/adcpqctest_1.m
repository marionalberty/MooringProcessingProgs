function [ib1,ib2,ib3,ib4,ib5,ib6] = adcpqctest_1(qcthresh,qc,u,w,erv,...
  moor_name,ADCPtype,sno,filen)

set(0,'defaultaxesfontsize',14,'defaultaxeslinewidth',0.7,...
  'defaultlinelinewidth',1,'defaultpatchlinewidth',0.7,...
  'defaultFigureColor','white')

% Modified by C. Germineaud (cyril.germineaud@legos.obs-mip.fr), Aug 2015
% Inputs: a structure of thresholds for each of the following:
%   qcthresh.errvel  :  error velocity
%   qcthresh.pgood   :  percent good from 4-beam solutions
%   qcthresh.cmag    :  correlation magnitude
%   qcthresh.vvel    :  vertical velocity
%   qcthresh.hvel    :  horizontal velocity
%   qcthresh.ea      :  echo amplitude

% Outputs: bin indexes for each test: ib1,...,ib6
% Called by: adcp_qc_comparison.m

err_vel = qcthresh.err_vel;  %test 1
pgood =qcthresh.pgood;   %test 2 Janet suggests 10 might be low
cmag = qcthresh.cmag;    %test 3
vvel = qcthresh.vvel;    % test 4
hvel = qcthresh.hvel;   %test 5
ea_thresh = qcthresh.ea_thresh;   %test 6

%notes
% 1. Error Velocity Test
% measurement of disagreement of measurement estimates of opposite beams.
% Derived from 2 idpt beams and therefore is 2 indp measures of vertical
% velocity
% pass if adcp.error_vel <= 15 cm/s
% suspect if 15 cm/s <error_vel < 30 cm/s
% fail if error_vel > 30 cm
% adcp.error_vel(bin,time) can have + or - values

% Test 2: percent good
% Percent good is the ratio of good pings for each ensemble.
% It ensures that the minimum number of samples and minimum theoretical
% standard deviation in the dataset are met.

% Test 3: correlation magnitude test
% pass if 3/4 correlation magnitude values for each bin are
% > 64 (WH LR 75kHz)
% > 110 (Ocean Observer 38 kHz NB)
% > 190 (Ocean Observer 38 kHz BB)
% in adcp.config.corr_threshold == 64
% but the values look bad when they are more than 110 so I'm going to go
% with that for now!

% test 4: vertical velocity test
% vertical velocity is average of 2 vertical velocities sensed by the 2
% beams, so test for 0.3 m/s threshold

% test 5: horizontal test,
% look for water speeds more than 1.25 m/s

% Test 6, Echo Amplitude test
% Echo Intensity: A key quality control parameter, echo intensity is a
% measure of the signal strength intensity returned to the transducer. High
% echo intensity can show solid targets (e.g. a boundary, obstruction or
% fish), while low echo intensity can show insufficient scatterers or the
% limits of profiling range for the environment.

% This test looks at the difference between consecutive bin values of ea
% and if the value exceeds the threshold, then the bin fails, and all bins
% above this are also considered to have failed.
% This test is only applied from the middle bin to the end bin, since it is
% a test designed to get rid of surface bins

%test 1, Error Velocity test
% measurement of disagreement of measurement estimates of opposite beams.
% Derived from 2 idpt beams and therefore is 2 indp measures of vertical
% velocity
% pass if adcp.error_vel <= 15 cm/s
% suspect if 15 cm/s <error_vel < 30 cm/s
% fail if error_vel > 30 cm

% adcp.error_vel(bin,time) can have + or - values and only one value
% given (not separate bins?)

ib1 = abs(erv) >= err_vel;

%test 2, Percent Good test for Long ranger, use only %good for 4 beam
% solutions (ie pg(4))
ib2 = qc(4).pg <= pgood;
% use 4 as it is the percentage of measurements that have 4 beam solutions

% Test 3, correlation magnitude test
% pass if 3/4 correlation magnitude values for each bin are
% > 64 (WH LR 75kHz)
% > 110 (Ocean Observer 38 kHz NB)
% > 190 (Ocean Observer 38 kHz BB)

% in adcp.config.corr_threshold == 64
% but the values look bad when they are more than 110 so I'm going to go
% with that for now!

isub1 = (qc(1).cr<=cmag);
isub2 = (qc(2).cr<=cmag);
isub3 = (qc(3).cr<=cmag);
isub4 = (qc(4).cr<=cmag);
% test nbins bins
isub_all = isub1+isub2+isub3+isub4;

% assign pass(0) or fail(1) values
% Where 3 or more beams fail, then the cmag test is failed
ib3 = isub_all > 2;

% Test 4, Vertical velocity test
ib4 = abs(w) >= vvel;

% Test 5, Horizontal velocity test
ib5 = abs(u) >= hvel;

%Find the number that fail the first five tests
ib = ib1 + ib2 + ib3 + ib4 + ib5;

[ifailar,~] = find(ib >= 2);

%Test 6, Echo Amplitude test
% this test looks at the difference between consecutive bin values of ea
% and if the value exceeds the threshold, then the bin fails, and all bins
% above this are also considered to have failed.
% This test is only applied from the middle bin to the end bin, since it is
% a test designed to get rid of surface bins

% note that csiro data is timexbins so need to invert this?
%
% [ii,jj] = size(u);
% ib6=zeros(size(u));
% ik = round(jj/2);
%
% ib6(:,ik+1:jj) = (diff(qc(1).ea(:,ik:jj),1,2)>ea_thresh)+ ...
%     (diff(qc(2).ea(:,ik:jj),1,2)>ea_thresh)+ ...
%     (diff(qc(3).ea(:,ik:jj),1,2)>ea_thresh)+ ...
%     (diff(qc(4).ea(:,ik:jj),1,2)>ea_thresh);

[ii,~] = size(u);
ib6=zeros(size(u));
ik = round(ii/2);

ib6(ik+1:ii,:) = (diff(qc(1).ea(ik:ii,:),1,1)>ea_thresh)+ ...
  (diff(qc(2).ea(ik:ii,:),1,1)>ea_thresh)+ ...
  (diff(qc(3).ea(ik:ii,:),1,1)>ea_thresh)+ ...
  (diff(qc(4).ea(ik:ii,:),1,1)>ea_thresh);

[ifailr,ifailc] = find(ib6 >= 1);
for ii = 1:length(ifailr);
  ib6(ifailr(ii):end,ifailc(ii))=1;
end
%ib6(ifailr:end,ifailc)=1;

[ifailr,~] = find(ib6 >= 1);
% diagnostic test
% inverted the ib as seems to plot all time series
% have commented this out while doing diagnostic plots

% % Tests 1 to 5: plot ib histogram
% hist(ib',[0,1,2,3,4,5])
% title(['Number of data that failed QC tests 1 to 5, by bin number - ' ...
%   moor_name '_' ADCPtype],'interpreter','None')
% 
% xlabel('Number of tests failed')
% colorbar

[ii,jj] = size(u);

totbin = 1:ii;    % from half bins to top bin

figure(1)
[N,~] = hist(ifailar,totbin);
% No is the number in the bin where other tests failed
barh(totbin,N);
% write in % in each bin where test contributed to other test failure
va = get(gca,'XLim');
x1 = va(end)/4;
for ibn = 1:length(totbin)
  num = ceil((N(totbin(ibn))/jj)*100);
  text(x1,totbin(ibn),[num2str(num) ' (' num2str(N(totbin(ibn))) ')'],...
    'color','r');
end
title(['Other Test Fails ' moor_name '_' ADCPtype '_' num2str(sno)],...
  'interpreter','None')
eval(['print -dpng ',[filen moor_name '_' ADCPtype '_' ...
  num2str(num2str(sno))],'_QCfailures1.png'])

% % Test 6: plot ib histogram
% hist(ib6',[0,1,2,3,4])
% title(['Number of data that failed QC test 6, surface bins - ' moor_name...
%   '_' ADCPtype],'interpreter','None')
% xlabel('Number of tests failed')
% colorbar

figure(2)
[N,~] = hist(ifailr,totbin);
% No is the number in the bin where other tests failed
barh(totbin,N);
% write in % in each bin where test contributed to other test failure
va = get(gca,'XLim');

x1 = va(end)/4;
for ibn = 1:length(totbin)
  num = ceil(N(totbin(ibn))/jj)*100;
  text(x1,totbin(ibn),[num2str(num) ' (' num2str(N(totbin(ibn))) ')'],...
    'color','r');
end

title(['Other EA Fails ' moor_name '_' ADCPtype '_' num2str(sno)],...
  'interpreter','None')
eval(['print -dpng ',[filen moor_name '_' ADCPtype '_' num2str(sno)],...
  '_QCfailures2.png'])

figure(11)
pcolor(ib)
shading flat
title(['Total QC Indicators: >= 2 is a bad bin ' moor_name '_' ADCPtype ...
  '_' num2str(sno)],'interpreter','None')
colorbar
eval(['print -dpng ',[filen moor_name '_' ADCPtype '_' num2str(sno)], ...
  '_totqcind.png'])

figure(12)
pcolor(ib6)
shading flat
colorbar
title(['Bins rejected by EA Criteria ' moor_name '_' ADCPtype '_' ...
  num2str(sno)],'interpreter','None')

eval(['print -dpng ',[filen moor_name '_' ADCPtype '_' num2str(sno)], ...
  '_totalqcea.png'])

% close specific figures
% careful, do not use clear all,
% figure(101) will be close too !
close(figure(1));close(figure(2));close(figure(11));close(figure(12))

end