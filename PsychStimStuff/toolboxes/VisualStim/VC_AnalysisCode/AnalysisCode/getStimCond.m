function [C StimCond ] = getStimCond(filename,bnoPrint)
% function [C TRIG_Cond] = getStimCond(filename,bnoPrint))
%
% Read in file with Stimulus Condition for each file
%
% INPUT: filename- cell containing 1 or several filenames
%        maxsample- sweeplength (all filenames must have same sweeplength)
%
% OUTPUT: C - list of all conditions
%         StimCond - condition for each Sweep/Trigger (Sweeps from consecutive files
%         are conconcatinated)
% Used with s_C_PSTH_units_multifile
%
% BA 071109
if nargin ==1; bnoPrint = 1; end
bYesAll = 1; 
nClast = []; StimCond = []; indNaN = [];

if ~iscell(filename)
    error('filename must be a cell')
end
for i = 1:length(filename)
%     if iscell(filename)
        stemp = filename{i}; 
%     else stemp = filename; end
    temp = load([stemp '_TrigCond']); % load trigger info
    % STIMFILE = load([PATH filename '_SFile']); % load stimulus filename info
    
    C = unique( temp.data(2,:)); % conditions
    C = C(find(~isnan(C))); % exclude NAN and deal with seperately
    
    if bnoPrint % report sweeps that have NaN Condition
        indNaN = find(isnan(temp.data(2,:))); printf('%s \t  %d of %d are NaN:',stemp, length(indNaN),size(temp.data,2)); 
        if ~isempty(indNaN); printf(' Sweeps: %s',num2str(indNaN)); 
        else printf(' Sweeps: none'); end
    end
    
    nC= length(C); % number of conditions
    
    % check all files have same number of conditions
%     if ~bYesAll
        if i>1 & nClast ~= nC
            stemp2 = ['datafiles have different number of Conditions. May be different stimulus:\n' num2str(nC) ' '  num2str(nClast)];
            warning(stemp2);
%             R = input([stemp2 '\n'  stemp ...
%                 '\n Continue? (Y/N/Y to (A)ll)'],'s');
%             switch upper(R)
%                 case 'N'
%                     error('terminated');
%                 case 'A'
%                     bYesAll = 1;
%             end
%         end
    end
    nClast = nC;
    
    StimCond = [StimCond temp.data(2,:)];
end


