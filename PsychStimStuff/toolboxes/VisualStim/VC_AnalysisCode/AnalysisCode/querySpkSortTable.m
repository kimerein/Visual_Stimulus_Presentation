function [statusout existFilenameSpkSort spkSortTable sexptname DAQchns] = querySpkSortTable(PATHEXPTLOG,sexptname,DAQchns,FilenameSpkSort,bupdate)
% function [statusout existFilenameSpkSort spkSortTable]  = querySpkSortTable(PATHEXPTLOG,sexptname,DAQchns,FilenameSpkSort,bAppend,bupdate)
% OR
% function [statusout existFilenameSpkSort spkSortTable]  = querySpkSortTable(PATHEXPTLOG,FilenameSpkSort)
% USAGE:
% for querying/updating kSortTable.
% When bupdate = 0 will return existFilenameSpkSort that matches sexptname
% and DAQchns specified.
%
% When bupdate = 1 and will empty entries that share sexptname and any of DAQchns
% then add an entry for sexptname with DAQchns provided
%
% IF 2 input arguments then this format
%   [statusout existFilenameSpkSort spkSortTable]  = querySpkSortTable(PATHEXPTLOG,FilenameSpkSort)
% ELSE
%   [statusout existFilenameSpkSort spkSortTable] = querySpkSortTable(PATHEXPTLOG,sexptname,DAQchns,FilenameSpkSort,bupdate)
%
% Table rows are unique by a combination of sexptname and DAQchns
% (no 2 entries should have the same sexptname and any single DAQchn value
% in common)
% INPUT:
% PATHEXPTLOG - path where exptTable is saved
% exptname string (e.g. concatenated DATE DEPTH MOUSE#)
% STF(N) struct with field 'filename'.
% associated with exptname.
% bupdate (optional default 0) if 1 will add FilenameSpkSort to sexptname
% OUTPUT
% statusout- NaN; % when bupdate = 1 status out is set to NaN
% statusout- -1  saved spkSortTable doesn't exit
% statusout- 0 exptname or FilenameSpkSort doesn't exit
% statusout- -2 exptname exists but sortfile with ALL daqchns does NOT exists
% statusout- 1 FOUND i.e. EITHER 1) sortfile for exptname with ALL daqchn exists (could include
%                                OR 2) FilenameSpkSort exists in table
% more daqchns not specified)
%
% existFilenameSpkSort - if statuout is 1 then this is the filename with
% spikesort data for the DAQchns
%
% spkSortTable is a N x 3  cell array
% #1 column is contains the exptname string
% #2 column has a cell contraining a number array of the DAQchns
% #3 column string with the Filename containing SpkSort data
% BA
% NOTE DOES NOT PROTECT YOU FROM PUTTING THE SAME spksortFilename IN MORE THAN ONE
% row
existFilenameSpkSort = [];
% if ~exist('bAppend','var'); bAppend=0; end % by default do not allow multiple entries with same exptname and DAQchn
if ~exist('bupdate','var'); bupdate=0; end

try
    filename = 'exptTable';
    load(fullfile(PATHEXPTLOG,filename),'spkSortTable');
    if ~exist('spkSortTable','var') % Table doesn't exist
        ind = 1;
        statusout = -1;
    else
        
        if nargin > 2
            % does expt exist already?
            ind = find(~cellfun(@isempty,strfind(spkSortTable(:,1),sexptname)));
            
            if ~isempty(ind) % if it does
                %         are are any of the same DAQchns used in sortfile with the same expt?
                % (Probably don't want to use DAQchns in more than one time in a
                % sortfile
                temp = cellfun(@cell2mat,spkSortTable(ind,2),'UniformOutput',false);
                temp = cellfun(@(x)ismember(DAQchns,x),temp,'UniformOutput',false);
                ind2 = find(cellfun(@(x)any(x==1),temp));
                if ~isempty(ind2) % expt entry exists where at least 1 DAQchn is used in sortfile
                    %                 if ~bAppend
                    spkSortTable(ind(ind2),:) = []; %  delete all sortentries that have the same DAQchns (if not append mode)
                    %                 end
                    
                    if length(ind2)~=length(DAQchns)
                        statusout = -2; % 1 or more DAQchns is not sorted in a file assocaite with exptname
                    else %all DAQchns are sorted in an file associated with exptname
                        statusout = 1; %
                        existFilenameSpkSort = spkSortTable{ind,3};
                    end
                    
                end
            else % exptname doesn't exist
                statusout = 0;
            end
            % find the first row that is entirely empty and put the new
            % Sort entry in there
            ind = find(cellfun(@isempty,spkSortTable),1,'first');
            if isempty(ind); ind = size(spkSortTable,1)+1; end
            
        elseif nargin ==2 % querySpkSortTable(PATHEXPTLOG,FilenameSpkSort)
            FilenameSpkSort = sexptname;
            ind = find(~cellfun(@isempty,strfind(spkSortTable(:,3),FilenameSpkSort)));
            if  ~isempty(ind) % if it does
                sexptname = spkSortTable{ind,1};
                DAQchns =  spkSortTable{ind,2};
                statusout = 1;
            else % not found
                statusout = 0;
                sexptname = NaN;
                DAQchns = NaN;
            end
        else
            error('wrong number of input arguments')
        end
        
        
    end
catch % create file if non exists
    ind = 1;
    filename = 'exptTable';
end



if bupdate
    try
        statusout = NaN;
        spkSortTable(ind,:) = {sexptname {DAQchns} FilenameSpkSort};
        save(fullfile(PATHEXPTLOG,filename),'spkSortTable','-append');
        
    catch ME
        getReport(ME)
        error('Must create exptTable before creating spkSortTable use queryFileExptTable')
    end
end