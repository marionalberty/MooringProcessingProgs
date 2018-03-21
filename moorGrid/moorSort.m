function DS = moorSort(DS)
% Sort the data within the intermediate data structures as a function of
% planned meter depth for each row.

% Get the indecies of the data sorted by planned meter depth
[~,iSS] = sort(DS.pdep);

% Get all the fields in the data structure
fnames = fieldnames(DS);

for i = 1:length(fnames)
  eval(sprintf('DS.%s = DS.%s(iSS,:);',fnames{i},fnames{i}))
end


%% If data structure includes velocity remove rows with no data

if sum(strcmp(fnames,'U')) > 0
  [~,t_len] = size(DS.U);
  percentNAN = sum(isnan(DS.U),2)/t_len;
  for i = 1:length(fnames)
    eval(sprintf('DS.%s(percentNAN > 0.95,:) = [];',fnames{i}))
  end
end
