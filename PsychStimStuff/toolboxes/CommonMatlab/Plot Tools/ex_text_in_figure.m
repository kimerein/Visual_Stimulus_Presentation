% example text in figure

subplot(stile(1),stile(2),4);
axis off
tmp = findstr(fn,'\');
stemp=sprintf('%s\n%s\n%s\n%s %s\n%d %d',date(),figdesc,fn(tmp(end)+1:end),strtrim(cell2mat(schan(1))),strtrim(cell2mat(schan(2))),bInh(1),bInh(2));
text(0,1,stemp,'Interpreter','None')
