function  [name] = getfname(filepath)
% function  [name] = getfname(filepath
tmp = findstr(filepath,'.');
tmp2 = findstr(filepath,'\');
name = filepath;
if ~isempty(tmp)
    tmp = max(tmp);
    if ~isempty(tmp2)
        name = filepath(max(tmp2)+1:tmp-1);
    else
        name = filepath(1:tmp-1);
    end
else
    if ~isempty(tmp2)
        name = filepath(max(tmp2)+1:end);
    end
end

