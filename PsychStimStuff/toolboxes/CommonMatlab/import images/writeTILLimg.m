function imdata = writeTILLimg(imdata,fname)
% function imdata = writeTILLimg(imdata,fname)
% saves F frames sequentially starting from the one specifed in fname
% in the tif format
% INPUT:
%        impdata R by C by F int16 of imagedata
%        fname: path and filename of first image file in sequence to write
%          fname must have the format <path><filename>_<#><.tif>
%extract number of first file
if ~isa(imdata, 'uint16')
    warning('imdata is NOT a uint16 (imdata will be type-cast)')
end
F = size(imdata,3);
initN = str2num(fname(max(strfind(fname,'_'))+1:max(strfind(fname,'.'))-1));
fend = fname(max(strfind(fname,'.')):end); % extension
fbeg= fname(1:max(strfind(fname,'_'))); %
for i=initN:initN+F-1
    cname = [fbeg num2str(i) fend];
    temp = double(imdata(:,:,i-initN+1));
%     save(cname,'temp','-ascii'); %ASCII save data (very large)
     imwrite(uint16(imdata(:,:,i-initN+1)),cname,'tif','Compression','none'); %load data
end