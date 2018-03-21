function s=sbe39readfunc(pr)
% sbe39pread.m  read an sbe39 temp  (possibly with pressure) recorder file
% return data in structure s
% % filenam=input('Enter filename: ','s');
% input file name as arg or get via window dialog
%
%  pr = deployment press.  Leave it off to read p from file
readpr=0;
if nargin == 0;
    readpr = 1;
end

[filenam, pathnam, filtindx]=uigetfile('*.asc','select sbe39 data input file ');
fid=fopen([pathnam filenam],'rt');

[fnam,fext]=strtok(filenam,'.');
diary([pathnam fnam '_rd.log']);
disp([datestr(now) ' Opening file ' filenam])
fid=fopen([pathnam filenam],'rt');
s.te=[];s.ti=[] ; s.pr=[];
if fid<=0
   disp(['could not open ' filenam])
   return
end

ii=1;
%temp=ones(10,1);
%time=ones(10,1);
data=0;
disp([datestr(now) ' reading and decoding file '   filenam])
tic;
hm=msgbox('starting...','Records Processed','none');drawnow
hm=msgbox('starting...','Records Processed','none','replace');
drawnow


nhdr = 5;
while 1
  line=fgetl(fid);
  nhdr = nhdr + 1;  
  if ~ischar(line)
     disp([datestr(now) ' *** '])
     fprintf('\nEnding at record number %d\n',ii);
     break
  end
  if strcmp(line,'*END*') 
       line=fgetl(fid);
        line=fgetl(fid);
         line=fgetl(fid);
          line=fgetl(fid);
         % skip 4 more lines  before good data
                    % # of header lines to skip when reading data
    break
  end
  
end

 fclose(fid);                        %  close file for subsequent textread of data
 
 

 if readpr
 
 [s.te,s.pr,da,mo,yr,hh,mm,ss]=textread([pathnam filenam],'%f ,%f ,%d %s %d,  %d:%d:%d%*[^\n]','headerlines',nhdr);
 
 
 else
  [s.te,da,mo,yr,hh,mm,ss]=textread([pathnam filenam],'%f ,%02d %s %4d,  %02d:%02d:%02d%*[^\n]','headerlines',nhdr);
      s.pr=pr*ones(size(s.te));
 end
colonsep=repmat(':',size(yr));
dashsep=repmat('-',size(yr));
spacesep=repmat(' ',size(yr));

          tim = [num2str(da,'%02d')  dashsep char(mo) dashsep num2str(yr) spacesep num2str(hh,'%02d') colonsep num2str(mm,'%02d')  colonsep num2str(ss,'%02d')];
         s.ti = datenum(tim);
    
			mess=['Done. ' num2str( length(yr) ) ' records.'];
        hm=msgbox(mess,'Records Processed','warn','replace');
        drawnow

t1=toc;
disp(['reading and decoding file took ' num2str(t1) ' secs'])
diary;

