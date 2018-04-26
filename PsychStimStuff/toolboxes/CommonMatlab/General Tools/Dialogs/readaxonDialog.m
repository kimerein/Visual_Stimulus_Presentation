function [fileheader filenumbers sSliceCellID] = readaxonDialog(bSliceCell)
%% function creates a dialog for user to enter filenheader and filenumbers
%% if bSliceCell = 1 then user will also be asked for the SliceCell
%% number
persistent answer;

strArray = java_array('java.lang.String',1);
strArray(1) = java.lang.String('Enter readfile HEADER:');
strArray(2) = java.lang.String('Enter readfile INDEX:');
num_lines = [1];
if bSliceCell
    strArray(3) = java.lang.String('Enter SliceCellID:');
%     num_lines = [3];
end
prompt = cell(strArray);
dlg_title = 'Read Filename'; 
if isempty(answer)
    defAns{1,1} = datestr(now, 'yyyy_mm_dd_');
    defAns{1,2} = '0';
    if (bSliceCell)
            defAns{1,3} = 'A1';
    end
else
    defAns = answer;
end
answer = inputdlg(prompt,dlg_title,num_lines,defAns);
fileheader = answer{1,1};

temp = strfind(answer{2,1},':')
if (~isempty(temp))
    filenumbers = [str2num(answer{2,1}(1:temp-1)):str2num(answer{2,1}(temp+1:end))];
else
    filenumbers = str2num(answer{2,1});
end

sSliceCellID = '-1';
if bSliceCell
   sSliceCellID = answer{3,1};
end
