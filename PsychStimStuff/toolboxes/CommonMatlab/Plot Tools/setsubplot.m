function setsubplot(fid,stile, spos,fhandle)
% function setsubplot(fid,stile, spos,fhandle)
% places figure with handle fhandle into a new figure (fid) in a subplot
% spos. 
%
% INPUT:
%      
%
av = getfig(fhandle);
if size(av) ~= size(spos)
    error('fhandle has more subplots than spos. i.e. There are more plots then locations specified')
end
figure(fid);
for i = 1:length(av)
    da = subplot(stile(1),stile(2),spos(i));
    na = copyobj(av(i),fid);
    set(na,'Position',get(da,'Position'));
    delete(da);
end