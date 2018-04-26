%% Load Sequence of pics (multiframe)
F = 10;% total number of frames
% R = 457; C = 535;
R = 688; C = 520;
% fname =  'G:\Scanziani Lab data (not Backup)\image data\s2020407nm7\Fluorescence 488nm_8.tif' ;% name of first file
% fname =  'G:\Scanziani Lab data (not Backup)\image data\s2020407nm9df\DeltaF Fluorescence 488nm_161.tif' ;% name of first file
fname =  'G:\Scanziani Lab data (not Backup)\image data\s2020407nm9\Fluorescence 488nm_10.tif' ;% name of first file
imdata9= loadTILLimg(fname,C,R,F);
%% subtract mean
mim = int16(mean(imdata9,3));
figure(1)
imagesc(mim)
colormap('gray')
im91 = double(imdata9(:,:,1)- mim)./double(mim);
% im91 = double(imdata9- repmat(mim,[1 1 size(imdata9,3)]))./mim;

%% df to compare
F = 1;% total number of frames
R = 688; C = 520;
fname =  'G:\Scanziani Lab data (not Backup)\image data\s2020407nm9df\DeltaF Fluorescence 488nm_10.tif' ;% name of first file
imf9= loadTILLimg(fname,C,R,F);

%% clunky ui
fid = 10;
cf =1;
imgui(imdata2,cf,fid)


imhist(imdata(:,:,1),2^16);

%% write images
% fname =  'G:\Scanziani Lab data (not Backup)\image data\s2020407nm7\df_8.tif' ;% name of first file
fname =  'G:\Scanziani Lab data (not Backup)\image data\s2020407nm7\df_8.txt' ;% name of first file
writeTILLimg(imdata2,fname);
