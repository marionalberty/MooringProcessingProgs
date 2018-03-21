function [qcthresh] = params_threshold_qc(params,beam_freq,aorient)
%====================================================================
% Key: SET UP THRESHOLD FOR ADCP QC TESTS

% Author: J. Sprintall (jsprintall@ucsd.edu), Dec 2007

% Edited and Modified as a function by: 
% C. Germineaud (cyril.germineaud@legos.obs-mip.fr), Aug 2015
%====================================================================

% Ouputs: a structure of thresholds for each of the following:
%   qcthresh.errvel  :  error velocity
%   qcthresh.pgood   :  percent good from 4-beam solutions
%   qcthresh.cmag    :  correlation magnitude
%   qcthresh.vvel    :  vertical velocity
%   qcthresh.hvel    :  horizontal velocity
%   qcthresh.ea      :  echo amplitude

if (beam_freq==300)
    % ADCP Workhorse 300 kHz parameters
    if strcmp(params.moor_name,'Solomon_M1') ||...
        strcmp(params.moor_name,'Solomon_M3') ||...
        strcmp(params.moor_name,'StGeorgesEast') ||...
        strcmp(params.moor_name,'StGeorgesWest')
        
        qcthresh.err_vel=0.15;  %test 1
        qcthresh.pgood=80;   %test 2
        qcthresh.cmag=110;   %test 3
        qcthresh.vvel=0.2;    % test 4
        qcthresh.hvel=2.0;   %test 5
        qcthresh.ea_thresh=30;   %test 6
    end
    
elseif (beam_freq==75)
    % ADCP Workhorse 75 kHz parameters
    if strcmp(params.moor_name,'Solomon_M1') ||...
        strcmp(params.moor_name,'Solomon_M2b') ||...
        strcmp(params.moor_name,'Solomon_M3') ||...
        strcmp(params.moor_name,'StGeorgesEast') ||...
        strcmp(params.moor_name,'StGeorgesWest') 
        
        qcthresh.err_vel=0.15;  %test 1
        qcthresh.pgood=50;   %test 2
        qcthresh.cmag=64; %test 3
        qcthresh.vvel=0.2;    % test 4
        qcthresh.hvel=2.0;   %test 5
        qcthresh.ea_thresh=25;   %test 6
        
    elseif strcmp(params.moor_name,'VitiazMiddle') &&...
            strcmp(aorient,'up') || strcmp(aorient,'down')
        qcthresh.err_vel=0.15;  %test 1
        qcthresh.pgood=50;   %test 2
        qcthresh.cmag=110;
        qcthresh.vvel=0.2;    % test 4
        qcthresh.hvel=2.0;   %test 5
        qcthresh.ea_thresh=25;   %test 6 
    end
else
    display('do not recognize this adcp type')
    pause
end

end

