function [spikesortdata rawdatainfo] = putStimCondSpkSortData(FilenameSpkSort)
% function [spikesortdata rawdatainfo swcond]  = putStimCondSpkSortData(FilenameSpkSort)
%
% function loads spikesortdata file (contains spiketimes and unit assigns for each spike)
% (generated with sortExtSpikes.m) and VisualStimCond for each sweep
% and as field to spikesortdata.file.VstimCond with the VisStimCond for eac spike
%
% Note even if file.VstimCond  exists will still reload the condition for
% each sweep and update the field
%
% BA
s_RIGSPECIFIC_SPIKESORT_CONFIG;
load(fullfile(SORTEDSPIKES_SAVEPATH,FilenameSpkSort),'spikesortdata','rawdatainfo');

% get fileheader of rawdata
for i= 1:length(rawdatainfo)
    [junk swcond] = getStimCond({fullfile(DAQSAVEPATH,rawdatainfo(i).filename)},0);
    temp= swcond( spikesortdata.file(i).sweep);
     if isrowvector(temp);         temp = temp';     end
     spikesortdata.file(i).VstimCond =temp;
end
    
save(fullfile(SORTEDSPIKES_SAVEPATH,FilenameSpkSort),'spikesortdata','-append')