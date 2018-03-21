% adcp_qc_comparison.m

% Programmer: Janet Sprintall
% Date: 06 December 2007


% Modified by Cyril Germineaud, for use on Solomon Sea moorings,
% put vars and params into structures and
% define params_threshold_qc.m function

% Date: August 2015

% Purpose: program compares various thresholds for the NDBC ADCP QC tests
% to determine what should be used for the ADCP data
% Idea is to retain more of the upper most data (don't want to throw out good surface
% data, but also don't want to retain the bad!)
% program will go through systematically and vary each parameter
% individually, which assumes each parameter acts individually. Which
% probably is not true for Error Velocity and % good, but probably is true
% for the other first 5 tests. The ea_thresh (test 6) is fairly indpt.

% 1. Want to know what % in each bin where the test value contributes to the
% failure of that cell
% 2. Want to know the % of each bin where the test value contributes to the
% failure of that cell, but not the ea (ib6)

% Assume standard set is by SIO ADCP data from NDBC
% adcp_file = [main_path filein moor_name SerNo '_sio.mat']; 
% main_path     main directory
% filein        subdirectory for that mooring
% moor_name     mooring name and deployment number
%               Dipolog_deploy1
% SerNo         serial number

% qcthresh =
%       err_vel: 0.1500
%         pgood: 10
%          cmag: 110
%          vvel: 0.3000
%          hvel: 2.000
%     ea_thresh: 30

close all;clear all;clc

set(0,'defaultaxesfontsize',14,'defaultaxeslinewidth',0.7,...
  'defaultlinelinewidth',1,'defaultpatchlinewidth',0.7,...
  'defaultFigureColor','white')

% Personnal paths:
driveName='/Users/marionsofiaalberty/MATLAB/Solomon_Sea/';
% driveName='/Users/cyrilgermineaud/Documents/MATLAB/';
% addpath([driveName 'Routines_Cyril'])
% addpath([driveName 'ToolBox'])

% Set up paths
params.main_path = [driveName 'Moorings/'];
% params.main_path = [driveName 'Solomon_Sea_Moorings/'];
params.file_dat = 'Data/';

% path for figures
params.figs_raw = [params.main_path 'Figures/Raw/QC_tests/'];

% Set up script options:
params.moor_depl=[1:8,10:12]; % mooring cases
params.print=1; % print figures
%%
% 1. INITIALIZE INPUT FOR EACH MOORING
%
% These parameters probably need to be changed for each mooring and
% deployment!
for idmoor=params.moor_depl
    disp 'on going mooring instrument...'
    disp(idmoor)
    
    % Clear variables
    clearvars -except driveName params idmoor
        
    % switch for mooring and deployment
    switch idmoor
        case 1
            
            % 1. St Georges East 300 khZ
            % --------------
            % for file names
            params.moor_name = 'StGeorgesEast';
            ADCPdir = 'ADCP/';
            % WH300 khz upward at 152 m
            
            WaterDepth = 1433;
            DeployDepth = 152;
            SerNo = 16832;
            ADCPtype = 'RDI300kHz';
            
            % adcp orientation
            aorient = 'up';
            lat = -4-6.174/60;
            lon = 152+31.116/60;
            
            % deployment times
            %start_time=datenum(2012,07,20,06,10,00);
            %end_time = datenum(2014,03,8,04,00,00);
            start_time=datenum(2012,07,20,07,04,24); % start data
            end_time = datenum(2014,03,8,04,34,24); % end data
            clock_drift = -11/3600; % convert into hrs
            
            % bin depths for plotting comparison purposes
            binz = 13:23;
            
        case 2
            
            % 1. St Georges East 75kHz
            % --------------
            
            % for file names
            params.moor_name = 'StGeorgesEast';
            ADCPdir = 'ADCP/';
            % LR75 khz downward at 154 m
            
            WaterDepth = 1433;
            DeployDepth = 154;
            SerNo = 16768;
            ADCPtype = 'RDI75kHz';
            
            % adcp orientation
            aorient = 'down';
            lat = -4-6.174/60;
            lon = 152+31.116/60;
            
            % deployment times
            %start_time=datenum(2012,07,20,06,10,00);
            %end_time = datenum(2014,03,8,04,00,00);
            start_time=datenum(2012,07,20,06,20,00); %start data
            end_time = datenum(2014,03,8,03,20,00); %end data
            clock_drift = 100/3600; % convert into hrs
            
            % bin depths for plotting comparison purposes
            binz = 25:35;
            
        case 3
            
            % 1. St Georges West 75kHz
            % --------------
            
            % for file names
            params.moor_name = 'StGeorgesWest';
            ADCPdir = 'ADCP/';
            
            WaterDepth = 1243;
            DeployDepth = 154;
            SerNo = 8866;
            ADCPtype = 'RDI75kHz';
            
            % adcp orientation
            aorient = 'down';
            lat = -4-6.174/60;
            lon = 152+33.804/60;
            
            % deployment times
            %start_time=datenum(2012,07,20,06,00,00);
            %end_time = datenum(2014,03,7,21,00,00);
            start_time=datenum(2012,07,20,02,11,12); % start data 
            end_time = datenum(2014,03,7,21,11,11);% end data
            clock_drift = -28/3600; % convert into hrs
            
            % bin depths for plotting comparison purposes
            binz = 25:35;
            
        case 4
            
            % 1. St Georges West 300kHz
            % --------------
            
            % for file names
            params.moor_name = 'StGeorgesWest';
            ADCPdir = 'ADCP/';
            
            WaterDepth = 1243;
            DeployDepth = 152;
            SerNo = 16833;
            ADCPtype = 'RDI300kHz';
            
            % adcp orientation
            aorient = 'up';
            lat = -4-6.174/60;
            lon = 152+33.804/60;
            
            % deployment times
            %start_time=datenum(2012,07,20,06,00,00);
            %end_time = datenum(2014,03,7,21,00,00);
            start_time=datenum(2012,07,20,01,32,24); % start data
            end_time = datenum(2014,03,7,20,32,24); % end data
            clock_drift = 69/3600; % convert into hrs
            
            % bin depths for plotting comparison purposes
            binz = 13:23;
            
        case 5
            
            % 5.  Vitiaz Middle 75kHz UP
            % --------------
            
            % for file names
            params.moor_name = 'VitiazMiddle';
            ADCPdir = 'ADCP/';
            
            % WH75 khz upward at 332 m
            WaterDepth = 1131;
            DeployDepth = 332;
            SerNo = 8998;
            ADCPtype = 'RDI75kHz';
            
            % adcp orientation
            aorient = 'up';
            lat = -5-56.65/60;
            lon = 147+46.68/60;
            
            % deployment times
            %start_time=datenum(2012,07,28,06,15,00);
            %end_time = datenum(2014,03,13,18,45,00);
            start_time=datenum(2012,07,28,23,00,00); % start data
            end_time = datenum(2014,02,14,10,59,59); % end data
            clock_drift = 375/3600; % convert into hrs
            
            % bin depths for plotting comparison purposes
            binz = 70:109;
            
        case 6
            
            % 6.  Vitiaz Middle 75kHz DOWN
            % --------------
            
            % for file names
            params.moor_name = 'VitiazMiddle';
            ADCPdir = 'ADCP/';
            
            WaterDepth = 1131;
            DeployDepth = 334;
            SerNo = 16811;
            ADCPtype = 'RDI75kHz';
            
            % adcp orientation
            aorient = 'down';
            lat = -5-56.65/60;
            lon = 147+46.68/60;
            
            % deployment times
            %start_time=datenum(2012,07,28,06,15,00);
            %end_time = datenum(2014,03,13,18,45,00);
            start_time=datenum(2012,07,29,00,10,24); % start data
            end_time = datenum(2014,03,14,01,10,24); % end data
            clock_drift = 124/3600; % convert into hrs
            
            % bin depths for plotting comparison purposes
            binz = 25:35;
            
        case 7
            
            %   Solomon Strait M1 300 kHz UP
            % ----------------
            
            % for file names
            params.moor_name = 'Solomon_M1';
            ADCPdir = 'ADCP/';
            
            WaterDepth = 2100;
            DeployDepth = 80;
            SerNo = 5307;
            ADCPtype = 'RDI300kHz';
            
            % adcp orientation
            aorient = 'up';
            lat = -4-57.48/60;
            lon = 153+6.024/60;
            
            % deployment times
            %start_time = datenum(2012,07,21,01,30,00);
            %end_time = datenum(2014,03,06,23,01,00);
            start_time = datenum(2012,07,20,23,24,00); % start data
            end_time = datenum(2014,03,06,19,24,00); % end data
            clock_drift = 540/3600; % convert into hrs
            
            % bin depths for plotting comparison purposes
            binz = 10:40;
            
        case 8
           
            % Solomon Strait M1 75 kHz DOWN
            % ----------------
            
            % for file names
            params.moor_name = 'Solomon_M1';
            ADCPdir = 'ADCP/';
            
            WaterDepth = 2100;
            DeployDepth = 102;
            SerNo = 3427;
            ADCPtype = 'RDI75kHz';
            
            % adcp orientation
            aorient = 'down';
            lat = -4-57.48/60;
            lon = 153+6.024/60;
            
            % deployment times
            %start_time = datenum(2012,07,21,01,30,00);
            %end_time = datenum(2014,03,06,23,01,00);
            start_time = datenum(2012,07,20,21,28,00); % start data
            end_time = datenum(2014,03,06,17,28,00); % end data
            clock_drift = 830/3600; % convert into hrs
            
            % bin depths for plotting comparison purposes
            binz = 1:30;
            
        case 9
            
            %         % Solomon Strait M2a 300 kHz UP
            %         % ----------------
            %
            %         % for file names
            %         params.moor_name = 'Solomon_M2a';
            %         ADCPdir = 'ADCP/';
            %
            %         WaterDepth = 2900;
            %         DeployDepth = 80;
            %         SerNo = 40005;
            %         ADCPtype = 'FQ300kHz'
            %
            %         % adcp orientation
            %         aorient = 'up';
            %         lat = -5-9.853/60;
            %         lon = 153+16.864/60;
            %
            %         % deployment times
            %         start_time = datenum(2012,07,16,04,17,00);
            %         end_time = datenum(2014,03,06,01,45,00);
            %
            %         % bin depths for plotting comparison purposes
            %         binz = 1:2;
            
        case 10
            
            % Solomon Strait M2b 75 kHz UP
            % ----------------
            
            % for file names
            params.moor_name = 'Solomon_M2b';
            ADCPdir = 'ADCP/';
            
            WaterDepth = 2710;
            DeployDepth = 400;
            SerNo = 1066;
            ADCPtype = 'RDI75kHz';
            
            % adcp orientation
            aorient = 'up';
            lat = -5-9.448/60;
            lon = 153+19.937/60;
            
            % deployment times
            %start_time = datenum(2012,07,15,22,11,00);
            %end_time = datenum(2014,03,06,02,24,00);
            start_time = datenum(2012,07,15,21,37,36); % start data
            end_time = datenum(2014,03,06,00,37,36); % end data
            clock_drift = 206/3600; % convert into hrs
            
            % bin depths for plotting comparison purposes
            binz = 1:30;
            
        case 11
            
            % Solomon Strait M3 300 kHz UP
            % ----------------
            
            % for file names
            params.moor_name = 'Solomon_M3';
            ADCPdir = 'ADCP/';
            
            WaterDepth = 2617;
            DeployDepth = 80;
            SerNo = 12143;
            ADCPtype = 'RDI300kHz';
            
            % adcp orientation
            aorient = 'up';
            lat = -5-8.507/60;
            lon = 154+17.975/60;
            
            % deployment times
            %start_time = datenum(2012,07,15,22,11,00);
            %end_time = datenum(2014,03,06,02,24,00);
            start_time = datenum(2012,07,17,00,50,24); % start data
            end_time = datenum(2014,03,04,15,50,24); % end data
            clock_drift = 774/3600; % convert into hrs
            
            % bin depths for plotting comparison purposes
            binz = 10:40;
            
        case 12
            
            % Solomon Strait M3 75 kHz DOWN
            % -----------------------------
            
            % for file names
            params.moor_name = 'Solomon_M3';
            ADCPdir = 'ADCP/';
            
            WaterDepth = 2617;
            DeployDepth = 102;
            SerNo = 14215;
            ADCPtype = 'RDI75kHz';
            
            % adcp orientation
            aorient = 'down';
            lat = -5-9.448/60;
            lon = 153+19.937/60;
            
            % deployment times
            %start_time = datenum(2012,07,15,22,11,00);
            %end_time = datenum(2014,03,06,02,24,00);
            start_time = datenum(2012,07,16,23,14,48); % start data
            end_time = datenum(2013,09,29,23,18,11); % end data
            clock_drift = 563/3600; % convert into hrs
            
            % bin depths for plotting comparison purposes
            binz = 20:30;
            
%         case 12
%             
%             % Solomon Strait M3 75 kHz DOWN part 1
%             % ----------------
%             
%             % for file names
%             params.moor_name = 'Solomon_M3';
%             ADCPdir = 'ADCP/';
%             
%             WaterDepth = 2617;
%             DeployDepth = 102;
%             SerNo = 14215.1;
%             ADCPtype = 'RDI75kHz';
%             
%             % adcp orientation
%             aorient = 'down';
%             lat = -5-9.448/60;
%             lon = 153+19.937/60;
%             
%             % deployment times
%             %start_time = datenum(2012,07,15,22,11,00);
%             %end_time = datenum(2014,03,06,02,24,00);
%             start_time = datenum(2012,07,16,23,14,48); % start data
%             end_time = datenum(2013,04,08,16,14,48); % end data
%             clock_drift = 563/3600; % convert into hrs
%             
%             % bin depths for plotting comparison purposes
%             binz = 20:30;
%             
%         case 13
%             
%             % Solomon Strait M3 75 kHz DOWN part 2
%             % ----------------
%             
%             % for file names
%             params.moor_name = 'Solomon_M3';
%             ADCPdir = 'ADCP/';
%             
%             WaterDepth = 2617;
%             DeployDepth = 102;
%             SerNo = 14215.2;
%             ADCPtype = 'RDI75kHz';
%             
%             aorient = 'down';
%             lat = -5-9.448/60;
%             lon = 153+19.937/60;
%             
%             % deployment times
%             %start_time = datenum(2012,07,15,22,11,00);
%             %end_time = datenum(2014,03,06,02,24,00);
%             start_time = datenum(2013,04,09,04,18,11); % start data
%             end_time = datenum(2013,09,29,23,18,11); % end data
%             clock_drift = 563/3600; % convert into hrs
%             
%             % bin depths for plotting comparison purposes
%             binz = 20:30;
    end
    
    %% 3. GET DATA INTO COMMON TIME BASE FOR START AND END DATE
    
    filein = ([params.main_path params.file_dat 'Raw/' params.moor_name '/' ADCPdir params.moor_name '_'...
        num2str(DeployDepth) 'm_' ADCPtype '_' num2str(SerNo) '.mat']);
    
    % Load filein
    load(filein)
    
    % Define ndate and depth variables
    time0=adcp.mtime;
    adcp.time = time0-clock_drift;
    ndate = adcp.time;
    depth = (adcp.depth' * ones(1,length(adcp.config.ranges)) - ones(length(adcp.depth),1)*adcp.config.ranges');
    
    % now get all times within the time period
    ii= ndate>=start_time&ndate<=end_time;
    depth=depth(ii,:);
    ndate=ndate(:,ii)';
    uu = adcp.east_vel(:,ii);
    vv = adcp.north_vel(:,ii);
    uorig = uu+i*vv;
    
    w = adcp.vert_vel(:,ii);
    erv = adcp.error_vel(:,ii);
    nbins = adcp.config.n_cells;

    % Raw data first check 
    plot_adcp_first_check(ii,params,ADCPtype,SerNo,DeployDepth,adcp)
    % close on-going figure
    close
    
    % now get qc parameters
    % set up a structure qc for diagnostics
    for nn =1:4
        qc(nn).ea = squeeze(adcp.intens(:,nn,ii));
        qc(nn).pg = squeeze(adcp.perc_good(:,nn,ii));
        qc(nn).cr = squeeze(adcp.corr(:,nn,ii));
    end
    
    %% 2. SET UP THRESHOLD FOR ADCP QC TESTS
    
    % SIO Standard Tests
    qcthresh.err_vel=0.15;  %test 1
    qcthresh.pgood=50;   %test 2 
    qcthresh.cmag=110;
    qcthresh.vvel=0.2;    % test 4
    qcthresh.hvel=2.0;   %test 5
    qcthresh.ea_thresh=30;   %test 6

    %% Do QC and velocity for qctest
    % for test switch
    for itt = 1:6
        test = itt;
        
        % limit to 4 different criteria so can plot easily
        
        switch test
            % CHANGE TEST HERE
            % test 1
            case 1
                qctest_type = 'err_vel';
                qctest_var = erv;
                qctest_val = [0.1 0.15 0.2 0.3];
                qctest = 1;
                
                % test 2
            case 2
                qctest_type = 'pgood';
                qctest_val = [10 30 50 80];
                qctest = 2;
                qctest_var = qc(4).pg;
                
                % test 3
            case 3
                % > 64 (WH LR 75kHz)
                % > 110 (Ocean Observer 38 kHz NB)
                % > 190 (Ocean Observer 38 kHz BB)
                % nanmean for WH300 is about 125
                qctest_type = 'cmag';
                qctest_val = [50 64 75 110];
                qctest = 3;
                qctest_var = qc(4).cr;
                
                % test 4
            case 4
                qctest_type = 'vvel';
                qctest_val = [0.15 0.2 0.3 0.35];
                qctest = 4;
                qctest_var = w;
                
                % test 5
            case 5
                qctest_type = 'hvel';
                qctest_val = [1.5 1.75 2.0 2.5];
                qctest = 5;
                qctest_var = abs(uorig);
                
                % test 6
            case 6
                qctest_type = 'ea_thresh';
                qctest_val = [20 25 30 35]; %[20 25 30 35];
                qctest = 6;
                qctest_var = (qc(1).ea);
        end
        
        
        % Set figure params
        figure(101)
        set(gcf,'PaperUnits','centimeters')
        
        for itest = 1:length(qctest_val);
            u = uorig;
            
            eval(['qcthresh.' qctest_type ' = qctest_val(itest)']);
            
            fprintf('Threshold Parameters for ADCP quality control\n\n')
            fprintf('Error Velocity %5.2f meters per second\n', qcthresh.err_vel)
            fprintf('Percent Good %5.2f \n',qcthresh.pgood)
            fprintf('Correlation Magnitude %5.2f \n',qcthresh.cmag)
            fprintf('Vertical Velocity %5.2f meters per second\n',qcthresh.vvel)
            fprintf('Speed %5.2f meters per second\n',qcthresh.hvel)
            fprintf('Echo Amplitude bin difference %5.2f \n\n',qcthresh.ea_thresh)
            
            % Define path for RAW figures per mooring
            params.subdir_raw = [params.figs_raw params.moor_name '/'];
            
            % 5. RUN DIAGNOSTICS FOR ADCP QC CHECK
            
            disp('calling adcp qc test')

            [ib1,ib2,ib3,ib4,ib5,ib6] = adcpqctest_1(qcthresh,qc,u,w,erv,...
                params.moor_name,ADCPtype,SerNo,params.subdir_raw);
            
            % put in big array so you can use it as a pivot test
            ib(1,:,:) = ib1;
            ib(2,:,:) = ib2;
            ib(3,:,:) = ib3;
            ib(4,:,:) = ib4;
            ib(5,:,:) = ib5;
            ib(6,:,:) = ib6;
            
            ib_others = ib1+ib2+ib3+ib4+ib5;
            
            % indices where cell failed other tests
            [ifailor,ifailoc] = find(ib_others>=2);
            
            theone = squeeze(ib(qctest,:,:));
            
            % these indicies contain those cells where the qctest case contributed to
            % the failure
            [ifailar,ifailac] = find((ib_others >= 2)&(theone==1));
            
            % failed ea test
            [ifailr,ifailc] = find(ib6 >= 1);
            
            for ii = 1:length(ifailr);
                ib6(ifailr(ii):end,ifailc(ii))=1;
            end
            
            [ifailr,ifailc] = find(ib6 >= 1);
            
            % now find the case where these cells caused the failure of the cell vs the
            % more pervasive ea test!
            [ifailnr,ifailnc] = find((ib_others >= 2)&(theone>=1)&(ib6==0));
            
            % ifailor = failed other tests
            % ifailar = failed other tests and test case contributed to that failure
            % ifailr  = failed ea test anyway
            % ifailnr = failed other tests and contributed, but did not fail ea test
            % ifailtr = failed cell total whether due to other tests or ea
            
            % Fig 1: want % in each bin where test contributed to other test failure
            % Fig 2: want % in each bin where test contributed to total failure but cell did not fail ea test
            
            % set bad data in velocity (ifailtr)
            % now total number of failed cells
            % first those cells that failed others test
            
            [ifailtr,ifailtc] = find(ib_others >= 2);
            ibr = ifailtr;  %
            ibc = ifailtc;
            
            for ir = 1:length(ibr)
                u(ibr(ir),ibc(ir)) = complex(nan,nan);
            end
            
            % now those cells that failed ea test
            ibr = ifailr;  %
            ibc = ifailc;
            
            for ir = 1:length(ibr)
                u(ibr(ir),ibc(ir)) = complex(nan,nan);
            end
            
            % now find what is the total number of failures is
            [ifailtr,ifailtc] = find(isnan(real(u)));
            
%             % For Lombok interested in plotting meridional velocity
%             if (params.moor_depl <= 2)
%                 is_v= imag(u);
%             else
%                 is_v = real(u);
%             end

            % Define is_v as the zonal velocity
            is_v = real(u);
            
            % histogram of number of rejected data in binz
            
            % ifailor = failed other tests
            % ifailar = failed other tests and test case contributed to that failure
            % ifailnr = failed other tests and contributed, but did not fail ea test
            % ifailr  = failed ea test anyway
            % ifailtr = total number of failed cells (the bottom line!)
            
            % Fig 1: want % in each bin where test contributed to other test failure
            % Fig 2: want % in each bin where test contributed to failure but cell did not fail ea test
            
            % Define totbin variable
            totbin = 1:size(u,1);
            
            if qctest~=6
                % only do for tests 1-5 which are the other tests!
                subplot(4,3,((itest-1)*3)+1)
                
                [No,xo] = hist(ifailor,totbin);     % No is the number in the bin where other tests failed
                [N,x] = hist(ifailar,totbin);    % N is number in the bin where test case contributed to other test failure
                barh(binz,N(binz));
                % write in % in each bin where test contributed to other test failure
                va = get(gca,'XLim');
                x1 = va(end)/4;
                for ibn = 1:length(binz)
                    num = ceil(N(binz(ibn))/No(binz(ibn))*100);
                    text(x1,binz(ibn),[num2str(num) ' (' num2str(N(binz(ibn))) ')'],'fontsize',10,'color','r');
                end
                title('% Test Case in OtherTest Fails') % QC: ' strrep(qctest_type,'_',' ') num2str(qctest_val(itest))])
                
                % ifailor = failed other tests
                % ifailar = failed other tests and test case contributed to that failure
                % ifailnr = failed other tests and contributed, but did not fail ea test
                % ifailr  = failed ea test anyway
                % ifailtr = failed cell through other and ea - the bottom line!
                
                % Fig 1: want % in each bin where test contributed to other test failure
                % Fig 2: want % in each bin where test contributed to failure but cell did not fail ea test
                
                subplot(4,3,((itest-1)*3)+2)
                [No,xo] = hist(ifailtr,totbin);     % No is the number in the bin where all test failed
                [N,x] = hist(ifailnr,totbin);    % N is number in the bin where test case contributed to other test failure and did not fail ea test
                barh(binz,N(binz));
                % write in % in each bin where test contributed to other test failure
                va = get(gca,'XLim');
                x1 = va(end)/4;
                for ibn = 1:length(binz)
                    num = ceil(N(binz(ibn))/No(binz(ibn))*100);
                    text(x1,binz(ibn),[num2str(num) ' (' num2str(N(binz(ibn))) ')'],'fontsize',10,'color','r');
                end
                
                set(gca,'tickdir','out')
                title('% Test Case (and not EA) in Total Fails') % QC: ' strrep(qctest_type,'_',' ') num2str(qctest_val(itest))])
            else     % case where test==6
                
                % test against the total number of cells
                subplot(4,3,((itest-1)*3)+2)
                [No,x] = hist(ifailtr,totbin);     % N is failed both ea and other tests
                [N,x] = hist(ifailr,totbin);    % N is number in the bin that failed ea test
                barh(binz,N(binz));
                % write in % in each bin where test contributed to other test failure
                va = get(gca,'XLim');
                x1 = va(end)/4;
                for ibn = 1:length(binz)
                    num = ceil(N(binz(ibn))/No(binz(ibn))*100);
                    text(x1,binz(ibn),[num2str(num) ' (' num2str(N(binz(ibn))) ')'],'fontsize',10,'color','r');
                end
                
                set(gca,'tickdir','out')
                title('Failed ea in Total Fails')
                
                
            end % other tests
            % CHANGE X,Y FOR EACH DEPLOYMENT
            
            if itest==1
                % title if itest = 1
                y1 = binz(end) + 10;
                %text(0,y1,[' Standard QC Test: ' params.moor_name],'fontsize',13)
            end
            
            % quality controlled velocity
            subplot(4,3,((itest-1)*3)+3);
            pcolor(ndate,binz,is_v(binz,:))
            shading flat
            hold on
            set(gca,'tickdir','out')
            
            % set title
            title(['  Test ' strrep(qctest_type,'_',' ') ' = ' num2str(qctest_val(itest))],'fontsize',13)
        end
        
        orient tall
        
        % save in the qc test directory
        nfop= sprintf('%s%s/%s_%s_%s_%s',params.figs_raw,params.moor_name,...
                      params.moor_name,ADCPtype,num2str(SerNo),qctest_type);
        nfoppdf=sprintf('%s.pdf',nfop);
        set(gcf,'units','centimeters')
        set(gcf,'papersize',[30 60])
        set(gcf,'paperposition',[0,0,30,60])
        print(gcf,'-dpdf','-r200',nfoppdf);
        
        % test for ea and then difference
        if qctest==6
            qctest_var = diff(qctest_var,1,1);
            totbin = totbin(1)+0.5:totbin(end)-0.5;
        end
        %% plot the mean value for that test parameter
        figure(102)
        clf
        set(gcf,'PaperUnits','centimeters')
        xSize = 20; ySize = 14;
        xLeft = (20-xSize)/2; yTop = (14-ySize)/2;
        set(gcf,'Position',[xLeft yTop xSize*40 ySize*40])
        set(gca,'fontsize',20)
        set(gca,'fontweight','bold')
        
        plot((nanmean(qctest_var')),totbin,'b','linewidth',3);
        hold on
        plot((nanmean(qctest_var'))+nanstd((qctest_var')),totbin,'r','linewidth',2);
        plot((nanmean(qctest_var'))-nanstd((qctest_var')),totbin,'r','linewidth',2);
        plot(max((qctest_var')),totbin,'g','linewidth',2);
        plot(min((qctest_var')),totbin,'g','linewidth',2);
        
        title([qctest_type ' Mean, Std, Min, Max: ' params.moor_name],...
            'interpreter','None','fontsize',15,'fontweight','bold')
        
        % draw on lines of thresholds
        for i = 1:length(qctest_val)
            plot([qctest_val(i) qctest_val(i)],[totbin(1) totbin(end)],'--c','linewidth',2);
        end
        
        % save in the qc test directory
        nfop= sprintf('%s%s/%s_%s_%s_%s_mean',params.figs_raw,params.moor_name,...
                      params.moor_name,ADCPtype,num2str(SerNo),qctest_type);
        nfoppng=sprintf('%s.png',nfop);
        print('-dpng','-r200',nfoppng);
        
        % save in the qc test directory
        %eval(['print -dpng ',[params.figs_raw params.moor_name '/' params.moor_name  '_' qctest_type],'_mean.png'])
        
    end % test switch
end