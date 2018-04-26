function imdata = loadTILLimg(fname,R,C,F)
% function loadTILLimg(fname,R,C,F)
% loads F frames sequentially starting from the one specifed in fname
% INPUT:
%        fname: path and filename of first image file in sequence to read
%                    fname must have the format <path><filename>_<#><.tif>
%        F: total number of frames
%        R: # rows
%        C: # col
% OUTPUT: imdata: R by C by F int16 of imagedata
%
% BA 020507

%extract number of first file
initN = str2num(fname(max(strfind(fname,'_'))+1:max(strfind(fname,'.'))-1));
fend = fname(max(strfind(fname,'.')):end); % extension
fbeg= fname(1:max(strfind(fname,'_'))); %

if nargin >=4
    imdata = zeros(R,C,F,'uint16');%predeclare
end
for i=initN:initN+F-1
    cname = [fbeg num2str(i) fend];
    imdata(:,:,i-initN+1) = imread(cname); %load data
end