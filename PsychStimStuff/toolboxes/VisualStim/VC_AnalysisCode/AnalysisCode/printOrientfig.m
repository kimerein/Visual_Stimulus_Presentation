function printOrientfig(bprintfile,fid,file)
% function printOrientfig(bprintfile,fid,file)
% bprintfile = 1 otherwise prints to printer
% use to print Orientation figure from analUnitTuningNEW.m
if nargin==0; bprintfile = 0; end
if nargin<=1; fid = gcf; end
if nargin <=2; file = 'Orientfig'; end

% file='largetmp2';

orient landscape;
P = get(fid,'position');
set(fid,'position',[0           0        3364*2         987*2])
% I do this stretching because subplot proportion of a figure is not
% contant with figure size (both on the display and printed)
% so makeing the figure big makes the subplots bigger when printed.
% this is independent of print or saveas commands and seems to be
% independent of what format figure is printed to
%
% for some reason there is a max on the size you can make a figure (this
% these numbers exceed the max on my 2 screen 1680 x1050 (don't know if it is
% resolution related

% dpi=300;
% print(gcf,'-dpng',['-r' num2str(dpi)],file);
if bprintfile
print(gcf,'-depsc2',file);
else
 print(gcf,'-dwinc');
end
   
set(fid,'position',P);

