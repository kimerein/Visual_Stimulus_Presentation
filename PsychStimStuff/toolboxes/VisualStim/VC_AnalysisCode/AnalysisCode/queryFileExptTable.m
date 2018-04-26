function   [statusout exptTable]= queryFileExptTable(PATHEXPTLOG,sexptname,STF,bupdate)
% function [statusout exptTable] = queryFileExptTable(PATHEXPTLOG,sexptname,STF,bupdate)
% for query ExptTable with  DAQfile names that are associate with an
% experiment.
% if bupdate is 1 then exptTable will be changed to included STF filenames
%
% INPUT:
% PATHEXPTLOG - path where exptTable is saved
% exptname string (e.g. concatenated DATE DEPTH MOUSE#)
% STF(N) struct with field 'filename'.
% OUTPUT
% exptTable is a N x 2  cell array
% #1 column is unique ID contains the exptname string
% #2 column has a cell contraining the strings of all DAQfilenames (from STF)
%     for this experiment
% statusout = -1 if exptname doesn't exist in table
% statusout = 1 if exptname with where STF are ALL of files in expt
% statusout = 2 if exptname with where STF are a SUBSUBSET of files
% statusout 0 i.e. 1 or more of STF are not in exptname
% BA
% NOTE DOES NOT PROTECT YOU FROM PUTTING THE SAME FILES IN MORE THAN ONE
% EXPT
   if nargin<4;  bupdate =1; end         

try
% %     temp =  dirc(fullfile(PATHEXPTLOG,'*.mat'),'ae','d');
% %     filename = temp{end,2}; % get most recent file
    filename = 'exptTable';
    load(fullfile(PATHEXPTLOG,filename),'exptTable');
    
    % does expt exist already?
%     ind = regexpcell(exptTable(:,1),sexptname);
    ind = find(~cellfun(@isempty,strfind(exptTable(:,1),sexptname)));

    if ~isempty(ind) % if it does
        % check if identical
        breplace = 0;
        if size(exptTable{ind,2},2)==size({STF.filename},2)
            temp = ismember({STF.filename},exptTable{ind,2});
            if any(temp==0) % any don't match
                breplace = 1 ;
            else
                statusout=1; % all values for exptname in table are in STF.filename
            end
        else
            temp = ismember({STF.filename},exptTable{ind,2});
            if all(temp)
                statusout = 2; % all STF.filename are within exptname (but expt has more files it)
            end
            breplace = 1;
        end
        if breplace % if not identical replace
            %     find first empty column or maxcolumn
            exptTable(ind,:) = [];
            exptTable(ind,:) = {sexptname {STF.filename}};
            disp('Existing Exptname replaced')
            % delete all entries in SpkSortTable with same exptname
            tableDeleteRowContains(fullfile(PATHEXPTLOG,filename),'spkSortTable',1,sexptname);
        else
            disp('Exptname with same filenames already exists')
            statusout = 1;
        end
    else % add expt
        statusout = -1; % exptname doesn't exit
        exptTable(end+1,:) = {sexptname {STF.filename}};
    end
catch % create file if non exists
    
    exptTable= {sexptname {STF.filename}};
    if bupdate
        save(fullfile(PATHEXPTLOG,filename),'exptTable');
    end
end
if bupdate
    save(fullfile(PATHEXPTLOG,filename),'exptTable','-append');
    statusout = 1;
end
