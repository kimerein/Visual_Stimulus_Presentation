function [statusout existFilenameSpkSort spkSortTable] = querySpkSortTable(PATHEXPTLOG,sexptname,DAQchns,FilenameSpkSort,bAppend,bupdate)
% function spkSortTable = querySpkSortTable(PATHEXPTLOG,sexptname,DAQchns,FilenameSpkSort,bAppend,bupdate)
% for querying/updating kSortTable.
% each row is unique by a combination of sexptname and DAQchns
% (no 2 entries should have the same sexptname and any single DAQchn value
% in common)
% INPUT:
% PATHEXPTLOG - path where exptTable is saved
% exptname string (e.g. concatenated DATE DEPTH MOUSE#)
% STF(N) struct with field 'filename'.
% CHECK THIS !!!!!!!!!!!! bAppend (optional default 0) if 1 Append filenames to current filenames
% associated with exptname.
% bupdate (optional default 0) if 1 will add FilenameSpkSort to sexptname
% OUTPUT
% statusout- NaN; % when bupdate = 1 status out is set to NaN
% statusout- -1  saved spkSortTable doesn't exit
% statusout- 0 exptname doesn't exit
% statusout- -2 exptname exists but sortfile with ALL daqchns does NOT exists
% statusout- 1 sortfile for exptname with ALL daqchn exists (could include
% more daqchns not specified)
%
% existFilenameSpkSort - if statuout is 1 then this is the filename with
% spikesort data for the DAQchns 
%
% exptTable is a N x 3  cell array
% #1 column is contains the exptname string
% #2 column has a cell contraining a number array of the DAQchns
% #3 column string with the Filename containing SpkSort data
% BA
% NOTE DOES NOT PROTECT YOU FROM PUTTING THE SAME Filename IN MORE THAN ONE
% row
existFilenameSpkSort = [];
if ~exist('bAppend','var'); bAppend=0; end
if ~exist('bupdate','var'); bupdate=0; end
try
    %     temp =  dirc(fullfile(PATHEXPTLOG,'*.mat'),'ae','d');
    %     filename = temp{end,2}; % get most recent file
    filename = 'exptTable';
    load(fullfile(PATHEXPTLOG,filename),'spkSortTable');
    if ~exist('spkSortTable','var') % Table doesn't exist
        ind = 1;
        statusout = -1;
    else
        %         ind = regexpcell(spkSortTable(:,1),sexptname);% does expt exist already?
        ind = find(~cellfun(@isempty,strfind(spkSortTable(:,1),sexptname)));
        
        if ~isempty(ind) % if it does
            %         are are any of the same DAQchns used in sortfile with the same expt?
            % (Probably don't want to use DAQchns in more than one time in a
            % sortfile
            temp = cellfun(@cell2mat,spkSortTable(ind,2),'UniformOutput',false);
            temp = cellfun(@(x)ismember(DAQchns,x),temp,'UniformOutput',false);
            
            ind2 = find(cellfun(@(x)any(x==1),temp));
            if ~isempty(ind2) % expt entry exists where at least 1 DAQchn is used in sortfile
                if ~bAppend % if not append delete all sortentries that have the same DAQchns
                    spkSortTable(ind(ind2),:) = [];
                end
                if length(ind2)~=length(DAQchns)
                    statusout = -2; % 1 or more DAQchns is not sorted in a file assocaite with exptname
                else %all DAQchns are sorted in an file associated with exptname
                    statusout = 1; %      
                    existFilenameSpkSort = spkSortTable(ind,3);
                end
            end
        else % exptname doesn't exist
             statusout = 0;
        end
        % find the first row that is entirely empty and put the new
        % Sort entry in there
        ind = find(cellfun(@isempty,spkSortTable),1,'first');
        if isempty(ind); ind = size(spkSortTable,1)+1; end
    end
catch % create file if non exists
    ind = 1;
    filename = 'exptTable';
end


spkSortTable(ind,:) = {sexptname {DAQchns} FilenameSpkSort};

if bupdate
try
     statusout = NaN; 
save(fullfile(PATHEXPTLOG,filename),'spkSortTable','-append');

catch ME
    getReport(ME)
    error('Must create exptTable before creating spkSortTable use queryFileExptTable')
end
end