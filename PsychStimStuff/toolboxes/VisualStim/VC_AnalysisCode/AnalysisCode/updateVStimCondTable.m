function [Vstimcellout VstimTable] = updateVStimCondTable(PATHEXPTLOG,filenameCell,bforceupdate)
% function [Vstimcellout VstimTable]  = updateVStimCondTable(PATHEXPTLOG,filenameCell,bforceupdate))
%
% load existing table of Visual parameters for a given DAQ filename
% this is a cell
%  - first row contains title information
%  - each subsequent rows contain value of visual stim parameters presented
% when a specific DAQ file was collected.
%  - First colum contains DAQnames
% which should be unique.
%
% %  e.g.
% 'DAQfilename'    'VarParm1'       'VarParmVal1'
% 'Data2009-11-09M1_012'    'Orientation'     [1x8 double]
%
% Input
%             filenameCellcell 1xN cell of filenames
%             PATHEXPTLOG  path to tables (usually defined in RIG specific files:TABLES_SAVEPATH
%
% Output
%             Vstimcellout - NxM  Cell with Vstim information for specified
%             filenames. N = length(filenameCell) + 1 (for the titles).
%             VstimTable - entire table of Vstims
% BA

if nargin <3; bforceupdate = 0; end % default don't force update, use Vstim from table if it exists.

DAQSAVEPATH = [];LOGFILEPATH = []; % predeclare because contains nested functions
S_RIGSPECIFIC_SPIKESORT_CONFIG;
filename = 'DAQfileConditionTables';

if ~iscell(filenameCell); error('filenames must be a cell'); end;


for i = 1:size(filenameCell,2) % for each daqfilename supplied
    bsave = 1 ; %internal flag
    if i==1 % load old Table
        try
            load(fullfile(PATHEXPTLOG,filename),'VstimTable');
            bnewtable = 0;
        catch % create file if non exists
            bnewtable = 1;
        end
    end


    % ** create cell array with parameters (first row titles, 2nd row data)
    clear tempcell;
    tempcell{1,1} = 'DAQfilename';
    tempcell{2,1} = filenameCell{i};

    if ~bnewtable || i>1
        % compare requested entry to current table by filename
        ind = find(cellfun(@(x) strcmp(x,tempcell(2,1)),VstimTable(:,1)));
        if ind & ~bforceupdate % if any match  and not updating
            tempcell = VstimTable([1 ind],:);  bsave = 0;
        elseif ind & bforceupdate % if match and UPDATING
            loadVParamhelper();
            if isequal(VstimTable(1,:),tempcell(1,:)) % check that columns names are the same
                %                 get rid of nans because nan~=nan
                tempcell1 = tempcell(2,:);
                tempcell1(cellfun(@any,cellfun(@isnan,tempcell1,'UniformOutput',false))) = {0};
                tempcell2 = VstimTable(ind,:);
                tempcell2(cellfun(@any,cellfun(@isnan,tempcell2,'UniformOutput',false))) = {0};
                temp = isequal(tempcell2,tempcell1);
                if temp % same entry exists
                    bsave = 0;
                else % entry is not same so replace
                    printf('***************************')
                    printf('DAQfile: %s with different parameters exists\n  ',tempcell{2,1});
                    printf('old entry:\n');
                    %                     VstimTable(ind,:)
                    printf('REPLACED with new:\n');
                    %                     tempcell(2,:)
                    VstimTable(ind,:) = tempcell(2,:);
                end
            else
                error('VstimTable Column Names do not match new entry column names')
            end
        else % doesn't exist already , load file and add new entry for this filename
            loadVParamhelper();
            VstimTable(end+1,:) = tempcell(2,:);
        end
    else % new table
        loadVParamhelper();
        VstimTable= tempcell;
    end


    if i ==1 %  output cell with Vstim info of filenames
        Vstimcellout = tempcell;
    else
        Vstimcellout(end+1,:) = tempcell(2,:);
    end
end


if bsave
    save(fullfile(PATHEXPTLOG,filename),'VstimTable');
end

function loadVParamhelper
                [VarParam Param] = getPsychStimParameters(filenameCell{i},DAQSAVEPATH,LOGFILEPATH);
                    % ** create cell array with parameters (first row
                    % titles, 2nd row data)
                for j = 1:2
                    if j>length(VarParam) % only 1 variable
                        tempcell{2,2+2*(j-1)} = NaN;
                        tempcell{1,2+2*(j-1)} = NaN;
                    else
                        tempcell{2,2+2*(j-1)} = VarParam(j).Name;
                        tempcell{2,3+2*(j-1)} = VarParam(j).Values;
                    end
                    tempcell{1,2+2*(j-1)}  =sprintf('VarParm%d',j);
                    tempcell{1,3+2*(j-1)}   =sprintf('VarParmVal%d',j);
                end
                
                temp = fieldnames(Param);
                indlocal = cellfun(@isempty, regexp(temp,'Additional'));% exclude unwanted fields
                
                % collect struct fieldsnames and values into first and 2nd row of cell
                tempval = struct2cell(Param);
                tempcell = [tempcell [temp(indlocal) tempval(indlocal)]'];
end
end