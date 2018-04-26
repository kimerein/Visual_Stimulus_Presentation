function addEdata(filepath,data)
% function addEdata(filepath,data);
% loads variable 'Edata' from filepath and sets
% Edata(i+1) = data
tmp = findstr(filepath,'\');

if (~isdir(filepath(1:tmp(end)))) %directory doesn't exist
    mkdir(filepath(1:tmp(end)));
end
if isempty(dir([filepath '.mat'])) % file doesn't exist
    i=0;
else
    load(filepath,'Edata')
    i = length(Edata);
end
if length(data) ==1
Edata(i+1) = data;
else 
    Edata = data
end
save(filepath,'Edata');