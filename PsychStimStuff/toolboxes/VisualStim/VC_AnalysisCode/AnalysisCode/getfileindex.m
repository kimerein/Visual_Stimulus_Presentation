function stemp = getfileindex(STF)
i = 1;
stemp = [];
while i <=length(STF); % get filename(s) and concatanate for savefile
    temp = findstr(STF(i).filename,'_');
    
    stemp = [stemp 'f' STF(i).filename(temp(end)+1:end)];% filename
    if i == 18  % impose a maximum length on filename
        i = length(STF)-2 ;
    end
    i = i+1;
end
end