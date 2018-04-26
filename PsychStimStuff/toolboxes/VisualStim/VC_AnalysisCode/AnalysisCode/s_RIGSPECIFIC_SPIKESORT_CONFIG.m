% all rig specific parameters for data analysis
% try and catch is a bit of a hack, but is necessary to make it possible to
% call this file from functions that have nested functions (where all
% variables must be defined when the function is initialized)
% to get around this one defines the desired parameter in the function e.g.
% DAQSAVEPATH = []; and then runs this file.
% all the undefined variables will throw errors which will be caught and
% the the defined ones will be set to the correct values...

try
DAQSAVEPATH = 'C:\Documents and Settings\Bass\My Documents\Scanziani Lab\VC proj\Data\Rawdata\';
catch
end
try
ANALYZED_DAQSAVEPATH = 'C:\Documents and Settings\Bass\My Documents\Scanziani Lab\VC proj\Data\Analyzed\';
catch
end
try
TABLES_SAVEPATH = 'C:\Documents and Settings\Bass\My Documents\Scanziani Lab\VC proj\Data\Analyzed\Tables\';
catch
end
try
SORTEDSPIKES_SAVEPATH = 'C:\Documents and Settings\Bass\My Documents\Scanziani Lab\VC proj\Data\Analyzed\SortedSpikes\';
catch
end

try
LOGFILEPATH = '\\132.239.203.21\My Documents\Matlab toolboxes\VStimLog\';
% LOGFILEPATH = 'G:\Bass_VstimComputer\VStimLog 112009\';
catch
end

try
VFIG_SAVEPATH = 'C:\Documents and Settings\Bass\My Documents\Scanziani Lab\VC proj\Data\Figures\';
catch
end

try
STA_SAVEPATH = 'C:\Documents and Settings\Bass\My Documents\Scanziani Lab\VC proj\Data\Analyzed\STA\';
catch
end

try
MOVIE_PATH = '\\132.239.203.21\My Documents\Matlab toolboxes\VStimConfig\movies\';
catch
end% global DAQSAVEPATH LOGFILEPATH TABLES_SAVEPATH ANALYZED_DAQSAVEPATH SORTEDSPIKES_SAVEPATH
