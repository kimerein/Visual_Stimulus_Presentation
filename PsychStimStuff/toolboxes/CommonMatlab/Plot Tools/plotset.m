function plotset(nsetup,fid)
% function plotset(nsetup,fid)
% setup plot according to config predefined nsetup
% dirc holding config files:

global DPATH
% if ~exist('dpath')
% %     dpath = '/Volumes/MatlabCode/DEFINE/';
%     dpath = 'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\MatLab\DEFINE\';
% end
sfn = ['c_plotsetup_' num2str(nsetup)];
if exist('fid')
    if ~isempty(fid)
        figure(fid);
    end
end
run([DPATH sfn]);
