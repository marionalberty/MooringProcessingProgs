% CdpSTR2Plot simple function to create the date (matlab time) and temperature from
% the .cdp file for Seabird 39s or 56s

% USAGE:
% >> CdpSTR2Plot
% Load the STR .cdp file using load and enter the STR filename into the

% OUTPUTS:
%     times and temps for plotting

% Oliver Vetter March 17th 2008
% updated: December 2012, Vetter
% changed the function to use textscan instead of fgetl loops.  much faster
% now.  also assigns the 4# serial number to the time and temp parameters

function [times,temp] = CdpSTR2Plot

[filename, pathname, ~] = uigetfile('*.cdp','Enter .cdp STR file'...
                                    ,'/Users/olivervetter/Documents/MATLAB/','multiselect','on');                                   
        
file = fullfile(pathname, filename);
fid = fopen(file);
sn = filename(end-16:end-13);

[param, ~] = textscan(fid,'%f%f%f%f%f%f%f','Delimiter',',');

times = datenum(param{1},param{2},param{3},param{4},param{5},param{6});
temp = param{7};

assignin('base',['time' sn],times);
assignin('base',['temp' sn],temp);
clear times temp

% 
% 
% %% 
% % Saving STR .mat files to text
% 
% %% 
% 
% %%%% OFU STR plot
% subplot 311
% plot(date1367,temp1367,'k', date3005,temp3005,'r');
%     datetick('x',12,'keepticks');set(gca,'xticklabel',[]);
%     l = legend('STR 1108','STR 1109',2);set(l,'box','off','orientation','horizontal');
%     axis([date3863(1) date3863(end) 26 32]);grid;
% subplot 312
% plot(date3544,temp3544,'k', date3036,temp3036,'r');
%     datetick('x',12,'keepticks');set(gca,'xticklabel',[]);ylabel('Temperature ^oC');
%     l = legend('STR 1152','STR 3034',2);set(l,'box','off','orientation','horizontal');
%     axis([date3863(1) date3863(end) 26 32]);grid;
% subplot 313
% plot(date3863,temp3863,'k', date1672,temp1672,'r');
%     datetick('x',12,'keepticks');xlabel('Date');
%     l = legend('STR 3038','STR 3039',2);set(l,'box','off','orientation','horizontal');
%     axis([date3863(1) date3863(end) 26 32]);grid; 
% 
% 
% %% TUT SST plot
% load ('V:\HA-10-01_Leg2 AmSamoa\DATA\Oceanography\Moorings\TUT\SST\TUT_SSTs.mat');
% subplot 211
% plot(dateSST268001,tempSST268001,'k');
%     set(gca,'xticklabel',[]);
%     l = legend('SST 268001',2);set(l,'box','off','orientation','horizontal');
%     axis([dateSST268002(1) dateSST268002(end) 26 32]);grid;
% subplot 212
% plot(dateSST268002,tempSST268002,'k');
%     datetick('x',12,'keepticks');xlabel('Date')
%     l = legend('SST 268002',2);set(l,'box','off','orientation','horizontal');
%     axis([dateSST268002(1) dateSST268002(end) 26 32]);grid;
%     ylabel('Temperature ^oC');
% 
% %% TUT STR plot
% subplot 611
% 
% plot(date1108,temp1108,'k', date1109,temp1109,'r');
%     datetick('x',12,'keepticks');set(gca,'xticklabel',[]);
%     l = legend('STR 1108','STR 1109',2);set(l,'box','off','orientation','horizontal');
%     axis([date1108(1) date1108(end) 26 32]);grid;
% subplot 612
% plot(date1152,temp1152,'k', date3034,temp3034,'r');
%     datetick('x',12,'keepticks');set(gca,'xticklabel',[]);
%     l = legend('STR 1152','STR 3034',2);set(l,'box','off','orientation','horizontal');
%     axis([date1108(1) date1108(end) 26 32]);grid;
% subplot 613
% plot(date3039,temp3039,'k', date3038,temp3038,'r');
%     datetick('x',12,'keepticks');set(gca,'xticklabel',[]);
%     l = legend('STR 3038','STR 3039',2);set(l,'box','off','orientation','horizontal');
%     axis([date1108(1) date1108(end) 26 32]);grid; ylabel('Temperature ^oC');
% subplot 614
% plot(date3538,temp3538,'k', date3541,temp3541,'r');
%     datetick('x',12,'keepticks');set(gca,'xticklabel',[]);
%     l = legend('STR 3538','STR 3541',2);set(l,'box','off','orientation','horizontal');
%     axis([date1108(1) date1108(end) 26 32]);grid;
% subplot 615
% plot(date3861,temp3861,'k', date3864,temp3864,'r');
%     datetick('x',12,'keepticks');set(gca,'xticklabel',[]);
%     l = legend('STR 3861','STR 3864',2);set(l,'box','off','orientation','horizontal');
%     axis([date1108(1) date1108(end) 26 32]);grid;
% subplot 616
% plot(date3865,temp3865,'k', date4605,temp4605,'r');
%     datetick('x',12,'keepticks');xlabel('Date')
%     l = legend('STR 3865','STR 4605',2);set(l,'box','off','orientation','horizontal');
%     axis([date1108(1) date1108(end) 26 32]);grid;
% 
% %% BAK SBE37 plot
% subplot 211
% plot(date,temp,'k');
%     datetick('x',12,'keepticks');set(gca,'xticklabel',[]);
%     l = legend('ODP 20m',2);set(l,'box','off');
%     axis([date(1) date(end) 24 32]);grid; ylabel('Temp ^oC');
%     title('Baker Subsurface Temperature Record');
% subplot 212
% plot(date,salinity,'k');
%     datetick('x',12,'keepticks');
%     axis([date(1) date(end) 33 36]);grid;
%     ylabel('Salinity PSU');xlabel('Date');
%         
% %% BAK STR plot
% subplot 411
% plot(date1649,temp1649,'k');
%     datetick('x',12,'keepticks');set(gca,'xticklabel',[]);
%     l = legend('STR 1649',2);set(l,'box','off');
%     axis([date1649(1) date1649(end) 24 32]);grid;
%     title('Baker Subsurface Temperature Record');
% subplot 412
% plot(date3081,temp3081,'k');
%     datetick('x',12,'keepticks');set(gca,'xticklabel',[]);
%     l = legend('STR 3081',2);set(l,'box','off');
%     axis([date1649(1) date1649(end) 24 32]);grid;
% subplot 413
% plot(date3088,temp3088,'k');
%     datetick('x',12,'keepticks');set(gca,'xticklabel',[]);
%     axis([date1649(1) date1649(end) 24 32]);grid;
%     l = legend('STR 3088',2);set(l,'box','off');
%     ylabel('Temperature');
% subplot 414
% plot(date4028,temp4028,'k');
%     datetick('x',12,'keepticks');
%     axis([date1649(1) date1649(end) 24 32]);grid;
%     l = legend('STR 4028',2);set(l,'box','off');
% 
% %% HOW STR plot
% subplot 411
% plot(date3085,temp3085,'k');
%     datetick('x',12,'keepticks');set(gca,'xticklabel',[]);
%     l = legend('STR 3085',2);set(l,'box','off');
%     axis([date3085(1) date3085(end) 24 32]);grid;
%     title('Howland Subsurface Temperature Record');
% subplot 412
% plot(date3087,temp3087,'k');
%     datetick('x',12,'keepticks');set(gca,'xticklabel',[]);
%     l = legend('STR 3087',2);set(l,'box','off');
%     axis([date3085(1) date3085(end) 24 32]);grid;
% subplot 413
% plot(date4024,temp4024,'k');
%     datetick('x',12,'keepticks');set(gca,'xticklabel',[]);
%     axis([date3085(1) date3085(end) 24 32]);grid;
%     l = legend('STR 4024',2);set(l,'box','off');
%     ylabel('Temperature');
% subplot 414
% plot(date4027,temp4027,'k');
%     datetick('x',12,'keepticks');
%     axis([date3085(1) date3085(end) 24 32]);grid;
%     l = legend('STR 4027',2);set(l,'box','off');
%    
% 
% %% JOH STR data plot
% subplot 511
% plot(date4025,temp4025,'k');
%     datetick('x',12,'keepticks');set(gca,'xticklabel',[]);
%     l = legend('STR 4025',2);set(l,'box','off');
%     axis([date4025(1) date4025(end) 23 30]);grid;
% subplot 512
% plot(date4026,temp4026,'k');
%     datetick('x',12,'keepticks');set(gca,'xticklabel',[]);
%     l = legend('STR 4026',2);set(l,'box','off');
%     axis([date4025(1) date4025(end) 23 30]);grid;
% subplot 513
% plot(date3914,temp3914,'k');
%     datetick('x',12,'keepticks');set(gca,'xticklabel',[]);
%     axis([date4025(1) date4025(end) 23 30]);grid;
%     l = legend('STR 3914',2);set(l,'box','off');
%     ylabel('Temperature');
% subplot 514
% plot(date4029,temp4029,'k');
%     datetick('x',12,'keepticks');set(gca,'xticklabel',[]);
%     axis([date4025(1) date4025(end) 23 30]);grid;
%     l = legend('STR 4029',2);set(l,'box','off');
% subplot 515
% plot(date3866,temp3866,'k');
%     datetick('x',12,'keepticks');%set(gca,'xticklabel',[]);
%     l = legend('STR 3866','Orientation','horizontal',2);set(l,'box','off');
%     axis([date4025(1) date4025(end) 23 30]);grid;
%     xlabel('Date'); 
% 
% %%
% plot(date,temp,'k');
%     datetick('x',12,'keepticks');
%     axis([date(1) date(end) 26.5 30.5]);
%     l = legend('Aunu''u Temp','Orientation','horizontal',2);set(l,'box','off');
%     grid;
%     ylabel('Temperature ^oC');xlabel('Date');
% 
% %% SAR
% plot(date1045,temp1045,'k',date1065,temp1065,'k');
%     datetick('x',12,'keepticks');
%     axis([date1065(1) date1045(end) 26 31]);
%     l = legend('SAR','Orientation','horizontal',1);set(l,'box','off');
%     grid;xx = get(gca,'xtick');
% ylabel('Temperature');
%      xlabel('Date'); 
%    
% 
% 
% %% TIN_AGU
% subplot 311
% plot(date3082,temp3082,'k');
%     datetick('x',12,'keepticks');set(gca,'xticklabel',[]);
%     axis([date3082(1) date3082(end) 27 31]);
%     l = legend('TIN 5.8 m','Orientation','horizontal',1);set(l,'box','off');
%     grid;xx = get(gca,'xtick');
% subplot 312
% plot(date1644,temp1644,'k');
%     datetick('x',12,'keepticks');set(gca,'xticklabel',[],'xtick',[xx]);
%     axis([date3082(1) date3082(end) 27 31]);
%     l = legend('TIN  13.4 m','Orientation','horizontal',1);set(l,'box','off');
%     grid;
%     ylabel('Temperature');
% subplot 313
% plot(date1564,temp1564,'k');
%     datetick('x',12,'keepticks');
%     axis([date3082(1) date3082(end) 27 31]);
%     set(gca,'xtick',[xx]);
%     l = legend('AGU 8.2m','Orientation','horizontal',1);set(l,'box','off');
%     grid;
%     xlabel('Date'); 
% 
% 
% %% LIS
% subplot 311
% plot(date1199,temp1199,'k',date3068,temp3068,'r');
%     datetick('x',12,'keepticks');set(gca,'xticklabel',[]);
%     axis([date1189(1) date1189(end) 18 32]);
%     l = legend('1199 0.5 m','3068 14 m','Orientation','horizontal',4);set(l,'box','off');
%     grid;
% subplot 312
% plot(date1189,temp1189,'k',date3073,temp3073,'r');
%     datetick('x',12,'keepticks');set(gca,'xticklabel',[]);
%     axis([date1189(1) date1189(end) 18 32]);
%     l = legend('1189 10 m','3073 20 m','Orientation','horizontal',4);set(l,'box','off');
%     grid;
%     xx = get(gca,'xtick');
%     ylabel('Temperature');
% subplot 313
% plot(date3079,temp3079,'k');
%     datetick('x',12,'keepticks');
%     axis([date1189(1) date1189(end) 18 32]);
%     set(gca,'xtick',[xx]);
%     l = legend('3079 23 m','Orientation','horizontal',4);set(l,'box','off');
%     grid;
%     xlabel('Date'); 
% 
% %% PHR
% 
% % South Shore Transect
% load PHR_STR_HI0809.mat
% subplot 311
% plot(date0164,temp0164,'k',date1048,temp1048,'r',date0362,temp0362,'b');
% datetick('x',12,'keepticks');set(gca,'xticklabel',[]);
% axis([date1048(1) date1048(end) 15 30]);
% l = legend('1048 1 m','0362 20m','0164 30m','Orientation','horizontal',4);set(l,'box','off');
% grid; ylabel('S Transect'); xx = get(gca,'xtick');
% subplot 312 
% plot(date0167,temp0167,'k',date0835,temp0835,'r',date0166,temp0166,'b');
%     datetick('x',12,'keepticks');set(gca,'xticklabel',[]);
%     axis([date1048(1) date1048(end) 15 30]);
%     l = legend('0167 22m','0835 22m','0166 23m','Orientation','horizontal',4);set(l,'box','off');
%     grid;
%    set(gca,'xtick',[xx]);
% ylabel('S Shore');
%     subplot 313
% plot(date0842,temp0842,'k',date1112,temp1112,'r');
%     datetick('x',12,'keepticks');
%     axis([date1048(1) date1048(end) 15 30]);
%     set(gca,'xtick',[xx]);
%     l = legend('0166 23m','1112 22m','Orientation','horizontal',4);set(l,'box','off');
%     grid;
%     xlabel('Date'); ylabel('S-E Shore');
% 
%  % NOrthshore and Lagoon data
%  % NOrthshore transect
%  subplot 311
% plot(date0361,temp0361,'k',date1195,temp1195,'r',date1196,temp1196,'b');
% datetick('x',12,'keepticks');set(gca,'xticklabel',[]);
% axis([date1048(1) date1048(end) 15 30]);
% l = legend('0361 1 m','1195 20m','1196 30m','Orientation','horizontal',4);set(l,'box','off');
% grid; ylabel('N Transect'); xx = get(gca,'xtick');
% subplot 312 
% plot(date0904,temp0904,'k',date1150,temp1150,'r');
%     datetick('x',12,'keepticks');set(gca,'xticklabel',[]);
%     axis([date1048(1) date1048(end) 15 30]);
%     l = legend('1368 22m','0904 22m','1150 23m','Orientation','horizontal',4);set(l,'box','off');
%     grid;
%    set(gca,'xtick',[xx]);
% ylabel('Lagoon');
%     subplot 313
% plot(date0835,temp0835,'k',date3246,temp3246,'r');
%     datetick('x',12,'keepticks');
%     axis([date1048(1) date1048(end) 15 30]);
%     set(gca,'xtick',[xx]);
%     l = legend('0835 1m','3246 1m','Orientation','horizontal',4);set(l,'box','off');
%     grid;
%     xlabel('Date');
% 
% %% KUR
% 
% % Plot KURE lagoon temps
% plot(date0161,temp0161,'k',date1149,temp1149,'k');
%     datetick('x',12,'keepticks');set(gca,'xticklabel',[]);
%     datetick
%     axis([date0161(1) date1149(end) 15 30]);
%     l = legend('Kure Lagoon N: 28.419, W: -178.345','Orientation','horizontal',4);set(l,'box','off');
%     grid;
%     xx = get(gca,'xtick');
% ylabel('Temperature');
% 
% load KUR_STR_SBE37.mat
% 
% subplot 311
% plot(date1049,temp1049,'k',date1149,temp1149,'r');
%     datetick('x',12,'keepticks');set(gca,'xticklabel',[]);
%     axis([date1149(1) date1149(end) 15 30]);
%     l = legend('1049 1 m','1149 10m','Orientation','horizontal',4);set(l,'box','off');
%     grid;
% subplot 312
% plot(date1190,temp1190,'k',date0836,temp0836,'r');
%     datetick('x',12,'keepticks');set(gca,'xticklabel',[]);
%     axis([date1149(1) date1149(end) 15 30]);
%     l = legend('1190 20 m','0836 1m','Orientation','horizontal',4);set(l,'box','off');
%     grid;
%     xx = get(gca,'xtick');
% ylabel('Temperature');
%     subplot 313
% plot(date37,temp37,'k',date0585,temp0585,'r',date37,salin37,'g');
%     datetick('x',12,'keepticks');
%     axis([date1149(1) date1149(end) 15   40]);
%     set(gca,'xtick',[xx]);
%     l = legend('SBE37','SBE39','Salinity','Orientation','horizontal',4);set(l,'box','off');
%     grid;
%     xlabel('Date'); ylabel('CREWS');
% 
% %% MID
% subplot 411
% plot(date1047,temp1047,'k');
%     datetick('x',12,'keepticks');set(gca,'xticklabel',[]);
%     axis([date1148(1) date1148(end) 15 30]);
%     l = legend('1047 1 m','Orientation','horizontal',4);set(l,'box','off');
%     grid;
% subplot 412
% plot(date1145,temp1145,'k');
%     datetick('x',12,'keepticks');set(gca,'xticklabel',[]);
%     axis([date1148(1) date1148(end) 15 30]);
%     l = legend('1145 3 m','Orientation','horizontal',4);set(l,'box','off');
%     grid;
%     xx = get(gca,'xtick');
% ylabel('Temperature');
%     subplot 413
% plot(date1148,temp1148,'k');
%     datetick('x',12,'keepticks');
%     axis([date1148(1) date1148(end) 15   30]);
%     set(gca,'xtick',[xx]);
%     l = legend('1148 1 m','Orientation','horizontal',4);set(l,'box','off');
%     grid;
% 
%     subplot 414
% plot(date37,temp37,'k',date37,salin37,'r');
%     datetick('x',12,'keepticks');
%     axis([date1148(1) date1148(end) 15   40]);
%     set(gca,'xtick',[xx]);
%     l = legend('Temp 1 m','Salinity','Orientation','horizontal',4);set(l,'box','off');
%     grid;
%     xlabel('Date'); ylabel('SBE 37');
% 
% %% LAY
% subplot 311
% plot(date1141,temp1141,'k');
%     datetick('x',12,'keepticks');set(gca,'xticklabel',[]);
%     axis([date1141(1) date1141(end) 20 30]);
%     l = legend('1141 1 m','Orientation','horizontal',4);set(l,'box','off');
%     grid;
% subplot 312
% plot(date1203,temp1203,'k');
%     datetick('x',12,'keepticks');set(gca,'xticklabel',[]);
%     axis([date1141(1) date1141(end) 20 30]);
%     l = legend('1203 3 m','Orientation','horizontal',4);set(l,'box','off');
%     grid;
%     xx = get(gca,'xtick');
% subplot 313
% plot(date1370,temp1370,'k');
%     datetick('x',12,'keepticks');
%     axis([date1141(1) date1141(end) 20 30]);
%     set(gca,'xtick',[xx]);
%     l = legend('1370 1 m','Orientation','horizontal',4);set(l,'box','off');
%     grid;
%     xlabel('Date');
% 
% %% MAR
% load MAR_HI0809_STRand37s
% subplot 411
% plot(date0834,temp0834,'k',date0583,temp0583,'r');
%     datetick('x',12,'keeplimits');set(gca,'xticklabel',[]);
%     axis([date0834(1) date0834(end) 20 30]);
%     l = legend('0834 CREWS Anchor 8m','0583 CREWS Buoy 0m','Orientation','horizontal',4);set(l,'box','off');
%     grid;
% subplot 412
% plot(date1142,temp1142,'k',date1143,temp1143,'r');
%     datetick('x',12,'keeplimits');set(gca,'xticklabel',[]);
%     axis([date0834(1) date0834(end) 20 30]);
%     l = legend('1142 4m','1143 1.5m','Orientation','horizontal',4);set(l,'box','off');
%     grid;
% subplot 413
% plot(date1146,temp1146,'k');
%     datetick('x',12,'keeplimits');set(gca,'xticklabel',[]);
%     axis([date0834(1) date0834(end) 20 30]);
%     l = legend('1146 14m','Orientation','horizontal',4);set(l,'box','off');
%     grid;
%    subplot 414
% plot(date37,temp37,'k',date37,salin,'r');
% datetick('x',12,'keeplimits');
%     axis([date0834(1) date0834(end) 20 37]);
%     l = legend('CREWS37 Temp 1.5m','CREWS37 Salinity 1m','Orientation','horizontal',4);set(l,'box','off');
%     grid;
%  xlabel('Date');
% 
% 
% %% FFS
% subplot 414
%     plot(date3036,temp3036,'k',date3084,temp3084,'r');
%     datetick('x',12,'keeplimits');
%     axis([date3036(1) date3036(end) 20 30]);
%     l = legend('3036 2m','Backreef North 3m','Orientation','horizontal');set(l,'box','off');
%     grid; xlabel('Date');
%     xx = get(gca,'xtick');
% subplot 411
%     plot(date1871,temp1871,'k',date1670,temp1670,'r');
%     datetick('x',12,'keeplimits');set(gca,'xticklabel',[],'xtick',[xx]);
%     axis([date3036(1) date3036(end) 20 30]);
%     l = legend('South Tran 25m','South Tran 10m','Orientation','horizontal');set(l,'box','off');
%     grid;
% subplot 412
%     plot(date3004,temp3004,'k');
%     datetick('x',12,'keeplimits');set(gca,'xticklabel',[],'xtick',[xx]);
%     axis([date3036(1) date3036(end) 20 30]);
%     l = legend('South Tran 2m','Orientation','horizontal');set(l,'box','off');
%     grid;ylabel('Temperature (C)');set(gca,'xticklabel',[]);
% subplot 413
%     plot(date1868,temp1868,'k',date3244,temp3244,'r');
%     datetick('x',12,'keeplimits');set(gca,'xticklabel',[],'xtick',[xx]);
%     axis([date3036(1) date3036(end) 20 30]);
%     l = legend('SST Anchor 10m','La Perouse','Orientation','horizontal');set(l,'box','off');
%     grid;
%   
% % 1670: middepth STR in south transect
% % 1868: SST Anchor NW Atoll
% % 1871: on EAR Anchor:  deep transect
% % 3004 : shallow 2m STR in transect
% % 3084: backreef 3m Whale-Skate island North
% % 3244: La Perouse Pinnacle