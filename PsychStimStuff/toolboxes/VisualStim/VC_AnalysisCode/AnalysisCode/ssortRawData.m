function [FilenameSpkSort] = ssortRawData(STF,exptparams,DAQchns,Extrparams,params)
% [FilenameSpkSort] =% ssortRawData(STF,exptparams,DAQchns,Extrparams,params)
% this function sorts STFs 
% this function is incomplete.. See TO DOs
%
% BA120109

s_RIGSPECIFIC_SPIKESORT_CONFIG;
sexptname = createExptname(STF,exptparams);

if ~exist('Extrparams','var')|| isempty(Extrparams)
    Extrparams.invt = 1;
Extrparams.TH_STD = 3;
end
if ~exist('params','var')|| isempty(params)
params.breSort = 0; % don't resort
params.bNoLoad = 0;% don't load previous extracted files
end

queryFileExptTable(TABLES_SAVEPATH,sexptname,STF,1);
[extractedSpkFilenames bReextracted] = extractDAQSpikes(STF,DAQchns,Extrparams,params.bNoLoad);

% MUST resort of reextracted 
if bReextracted; params.breSort = 1; end
% MUST FIX FILENAME saved in sortExSpike
% speed up plotting ISI and autocorrelation
[FilenameSpkSort reSorted] = sortExtSpikes(extractedSpkFilenames,params.breSort);
 % TO DO Fig durectory
 % coudl be faster to pass the spikesort struct and not have to load it
 % from disk again...
 % *** TO DO CARRY through deletes in side this function to UNIT TABLE
 if reSorted % only update the table if resorted    
     querySpkSortTable(TABLES_SAVEPATH,sexptname,DAQchns,FilenameSpkSort,1);
 end


