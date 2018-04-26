function [exptTable statusout] = queryFileExptTable(PATHEXPTLOG,sexptname,STF,bupdate)
% function exptTable = queryFileExptTable(PATHEXPTLOG,sexptname,STF,bupdate)
% for updateing ExptTable with new DAQfile names that are associate with an
% experiment.
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
%
% statusout = 1 if exptname with where STF are ALL of files in expt
% statusout = 2 if exptname with where STF are a SUBSUBSET of files
% otherwise 0 
% BA
% NOTE DOES NOT PROTECT YOU FROM PUTTING THE SAME FILES IN MORE THAN ONE
% EXPT
            

try
% %     temp =  dirc(fullfile(PATHEXPTLOG,'*.mat'),'ae','d');
% %     filename = temp{end,2}; % get most recent file
    filename = 'exptTable';
    load(fullfile(PATHEXPTLOG,filename),'exptTable');
    
    % does expt exist already?
    ind = regexpcell(exptTable(:,1),sexptname);
    if ~isempty(ind) % if it does
        % check if identical
        breplace = 0;
        if size(exptTable(ind,2),2)==size({STF.filename},2)
            if any(strcmp(exptTable{ind,2},{STF.filename})==0);
                breplace = 1 ;
            else
                statusout=2;
            end
        else
            temp = strcmp(exptTable{ind,2},{STF.filename});
            breplace = 1; end;
        
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
end
