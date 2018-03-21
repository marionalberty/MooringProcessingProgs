function logAns = SNexists(data,TT)
% Check if the intstrument in question has already been added to the
% intermediate temperature grid. Return logical true if it already exsists
% within the serial number cell array and logical false otherwise.

% Get instrument SN and used SNs
instSN = data.serialNo;
moorSN = TT.snum;

% Compare serial numbers
icomp = strcmp(moorSN,instSN);

if sum(icomp) > 0
  logAns = true;
else
  logAns = false;
end