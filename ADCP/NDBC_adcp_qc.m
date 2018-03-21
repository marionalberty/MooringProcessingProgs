function [ifailar, ifailac,ifailr,ifailc] = adcpqctest(qcthresh,qc,u,w,erv,moor_name,filen)

% Inputs: a structure of thresholds for each of the following:
%   qcthresh.errvel  :  error velocity
%   qcthresh.pgood   :  percent good from 4-beam solutions
%   qcthresh.cmag    :  correlation magnitude
%   qcthresh.vvel    :  vertical velocity
%   qcthresh.hvel    :  horizontal velocity
%   qcthresh.ea      :  echo amplitude

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

% test 2: percent good

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
% Echo Intensity: A key quality control parameter, echo intensity is a measure of the signal
% strength intensity returned to the transducer. High echo intensity can show solid targets
% (e.g. a boundary, obstruction or fish), while low echo intensity can show insufficient
%     scatterers or the limits of profiling range for the environment.

% This test looks at the difference between consecutive bin values of ea and
% if the value exceeds the threshold, then the bin fails, and all bins
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

%test 2, Percent Good test for Long ranger, use only %good for 4 beam solutions (ie pg(4))
ib2 = qc(4).pg <= pgood;  %use 4 as it is the percentage of measurements that have 4 beam solutions

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
[ifailar,ifailac] = find(ib >= 2);

%Test 6, Echo Amplitude test
% this test looks at the difference between consecutive bin values of ea and
% if the value exceeds the threshold, then the bin fails, and all bins
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

[ii,jj] = size(u);
ib6=zeros(size(u));
ik = round(ii/2);

ib6(ik+1:ii,:) = (diff(qc(1).ea(ik:ii,:),1,1)>ea_thresh)+ ...
    (diff(qc(2).ea(ik:ii,:),1,1)>ea_thresh)+ ...
    (diff(qc(3).ea(ik:ii,:),1,1)>ea_thresh)+ ...
    (diff(qc(4).ea(ik:ii,:),1,1)>ea_thresh);

[ifailr,ifailc] = find(ib6 >= 1);
ib6(ifailr:end,ifailc)=1;

[ifailr,ifailc] = find(ib6 >= 1);

% save adcpqc_ib.mat ib*

% diagnostic test
% inverted the ib as seems to plot all time series

% have commented this out while doing diagnostic plots
% % hist(ib,[0,1,2,3,4,5])
% % title(['Number of data that failed QC tests 1 to 5, by bin number - ' moor_name])
% % xlabel('Number of tests failed')
% % colorbar
% % eval(['print -djpeg ',[filen moor_name],'QCfailures1.jpg'])
% % 
% % hist(ib6,[0,1,2,3,4])
% % title(['Number of data that failed QC test 6, surface bins - ' moor_name])
% % xlabel('Number of tests failed')
% % colorbar
% % eval(['print -djpeg ',[filen moor_name],'QCfailures2.jpg'])

% figure(1)
% pcolor(ib)
% shading flat
% title('Total QC Indicators: >= 2 is a bad bin')
% colorbar
% eval(['print -djpeg ',[filen moor_name], '_adcp_totalqcind.jpg'])
% 
% figure(2)
% pcolor(ib6)
% shading flat
% colorbar
% title(['Bins rejected by EA Criteria Difference ',num2str(qcthresh.ea_thresh)])
% eval(['print -djpeg ',[filen moor_name], '_adcp_totalqcea.jpg'])


end