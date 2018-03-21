% read and sort the Temperature Logger data into sensible times and QC a
% little

close all
clear all

% common paths

addpath('/Applications/MATLAB74/toolbox/plots')

main_path = '/Users/jsprintall/Documents/JanetsDocuments/swpacific/MoorSPICE/Moorings/';

moor_depl = 1;

% path for figures
file_fig = 'Figures/'
file_dat = 'Data/'


%%
% 1. INITIALIZE INPUT FOR EACH MOORING
%
% These parameters probably need to be changed for each mooring and
% deployment!

switch moor_depl
case 1

        % ST Georges East
        
        moor_name = 'StGeorgesEast';
        mn = 'SGE'

        inst = 'SIO Temperature Logger';

         lat = -4-6.174/60;
        lon = 152+31.116/60;
        waterz = 1433;

        % deployment times
        start_time=datenum(2012,07,20,06,10,00);
        end_time = datenum(2014,03,8,04,00,00);

        sno{1} = '165';
        plannedz{1} = 130;
        fname{1} = 'SGE-0130m_SIO-TL165_Spice2012_19Apr2012.csa'
        
        sno{2} = '160'
        plannedz{2} = 180;
        fname{2} = 'SGE-0350m_SIO-TL147_Spice2012_19Apr2012.csa';
        
        sno{3} = '161';
        plannedz{3} = 1100;
        fname{3} = 'SGE-1100m_SIO-TL161_Spice2012_19Apr2012-01.csa';
        
        
        sno{4} = '161';
        plannedz{4} = 1100;
        fname{4} = 'SGE-1100m_SIO-TL161_Spice2012_19Apr2012-02.csa';
        

    case 2

        % ST Georges West
        
        moor_name = 'StGeorgesWest';
        mn = 'SGW'

        inst = 'SIO Temperature Logger';

        
        lat = -4-6.174/60;
        lon = 152+33.804/60;
        waterz = 1433;

        % deployment times
        start_time=datenum(2012,07,20,06,00,00);
        end_time = datenum(2014,03,7,21,00,00);
        
        sno{1} = '162';
        plannedz{1} = 80;
        fname{1} = 'SGW-0080m_SIO-TL162_Spice2012_19Apr2012.csa';
        
        sno{2} = '138';
        plannedz{2} = 180;
        fname{2} = 'SGW-0180m_SIO-TL138_Spice2012_19Apr2012.csa';
        
        sno{3} = '139';
        plannedz{3} = 800;
        fname{3} = 'SGW-0800m_SIO-TL139_Spice2012_19Apr2012.csa';


% % % % % %         %%
% % % % % %         % Deployment2 Dipolog
% % % % % %         %-----------------
% % % % % %         % for file names
% % % % % %         moor_name = 'Dipolog_deploy2_';
% % % % % %         % path for prints and other output in process subdirectory
% % % % % %         filen = 'deploy2/dipo/process/'
% % % % % %         % path for input file from raw subdirectory
% % % % % %         filenr = 'deploy2/dipo/raw/'
% % % % % %         fname = 'Dipolog_TL161_485m.csa'



   

end


%% load csa file
ninst = length(sno);

for nf = 1:ninst

filein = [main_path file_dat moor_name '/SIOT/' fname{nf}]

t = loadcsa(filein);

rmfield(t,'depth');
rmfield(t,'lat');
rmfield(t,'lon');


% save data within start and end times
isub = (t.times>=start_time&t.times<=end_time);
s = struct('times', t.times(isub));
s = setfield(s,'data', t.data(isub));
s= setfield(s,'startdate_mat',t.times(1));
s= setfield(s,'stopdate_mat', t.times(end));
s= setfield(s,'latitude', lat);
s= setfield(s,'longitude', lon);
s= setfield(s,'plannedMeterDepth', plannedz{nf});
s= setfield(s,'waterDepth', waterz);
s= setfield(s,'meterType',inst);
s= setfield(s,'serialNo', sno);
% % 




%% plot raw data
figure(1)
clf
set(gcf,'defaultaxesFontWeight','bold')
set(gcf,'defaultAxesFontSize',15)
set(gcf,'DefaultAxesLineWidth',2)

plot(s.times,s.data)

axis([s.times(1) s.times(end) min(s.data) max(s.data)])
montick('x','m',get(gca,'position'))

ylabel('Temperature')
title(['SIO TL' char(sno(nf)) ' ' moor_name ' z= ' num2str(s.plannedMeterDepth)])

outfile = [main_path file_fig moor_name '/' moor_name '_' (num2str(plannedz{nf})) 'm_' char(sno{nf}) '_TL.pdf']
eval(['print -dpdf ' outfile]);
%% save data

tl_file = [main_path file_dat moor_name '/SIOT/'  moor_name '_' num2str(plannedz{nf}) 'm_' char(sno{nf}) '_qc.mat'];

clear t
t=s;
eval(['save ' tl_file ' t']);

end