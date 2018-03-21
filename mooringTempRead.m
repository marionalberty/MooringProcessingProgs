function t=mooringTempRead(fname)
% Input:
% fname = a full filename for any raw temperature data file from the
% moorspice moorings
%
% Output:
% t = structure containing a header character array and a data matrix

% There are 4 file types possible for mooring temperature data: cnv, asc,
% csv, dat, and csa. This function takes the given file for a temperature
% time series and produces a numerical matrix of the timeseries data and
% a charater array with the file header.


% Open file for reading as a text file
fid=fopen(fname,'rt');
% Get file type
fext=fname(end-2:end);

% Initialize header cell
header=cell(1);
data=[];

i=1;
% Check which file type
if strcmp(fext,'cnv')
    % .cnv files are used for SBE-37 on French moorings
    
    % Make cell array of header strings
    str=fgetl(fid);
    % Stop when you reach end mark
    while ~strcmp(str,'*END*')
        if strncmp(str,'# nquan',7)
            nvars=str2double(str(end-1:end));
        end
    header{i,1}=str;
    i=i+1;
    str=fgetl(fid);
    end
    header{i,1}=str;
    % Make data matrix
    data=transpose(fscanf(fid,'%f',[nvars inf]));
    

elseif strcmp(fext,'asc')
    % .asc files are used for SBE-37 on SIO moorings and all SBE-39
    
    % Make cell array of header strings
    str=fgetl(fid);
    % Stop when you reach end mark
    while ~strncmp(str,'start sample number',19)
        header{i,1}=str;
        i=i+1;
        str=fgetl(fid);
    end
    header{i,1}=str;
    % Make data matrix
    dataR=textscan(fid,'%s','Delimiter',',');
    nvars=find(strncmp('00:',dataR{1},3),1);
    dataR=reshape(dataR{1},nvars,[])';
    date=char(dataR{:,nvars-1});
    time=char(dataR{:,nvars});
    dnum=datenum([date time],'dd mmm yyyyHH:MM:SS');
    for i=1:nvars-2
        data=[data str2num(char(dataR{:,i}))];
    end
    data=[data dnum];
    

elseif strcmp(fext,'csv')
    % .csv files are used for all SBE-56
    
    % Make cell array of header strings
    str=fgetl(fid);
    % Stop when you reach end mark
    while ~strncmp(str,'"Date"',6) && ~strncmp(str,'"Sample',7)
        header{i,1}=str;
        i=i+1;
        str=fgetl(fid);
    end
    header{i,1}=str;
    % Make data matrix
    if strncmp(str,'"Date"',6)
        dataR=textscan(fid,'%q %q %q','Delimiter',',');
        date=char(dataR{1});
        time=char(dataR{2});
        dnum=datenum([date time],'yyyy-mm-ddHH:MM:SS');
        data=[str2num(char(dataR{3})) dnum];
    else
        dataR=textscan(fid,'%q %q %q %q','Delimiter',',');
        date=char(dataR{2});
        time=char(dataR{3});
        dnum=datenum([date time],'yyyy-mm-ddHH:MM:SS');
        data=[str2num(char(strrep(dataR{4},',','.'))) dnum];
    end
    
    
elseif strcmp(fext,'csa')
    % .csa files are used for all SIOT
    
    % Make cell array of header strings
    str=fgetl(fid);
    % Stop when you reach end mark
    while ~strncmp(str,'deg_C',5)
        header{i,1}=str;
        i=i+1;
        str=fgetl(fid);
    end
    header{i,1}=str;
    % Make data matrix
    data=textscan(fid,'%f');
    data=data{1};
    
    
elseif strcmp(fext,'dat')
    % .dat files are used for all RBR instruments
    
    % Make cell array of header strings
    str=fgetl(fid);
    % Stop when you reach end mark
    while ~strcmp(strtrim(str),'Temp')
        header{i,1}=str;
        i=i+1;
        str=fgetl(fid);
    end
    header{i,1}=str;
    % Make data matrix
    data=textscan(fid,'%f');
    data=data{1};
    
    
end

fclose(fid);

t.data=data;
t.header=char(header);