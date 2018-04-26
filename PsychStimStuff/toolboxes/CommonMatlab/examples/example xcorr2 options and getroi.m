%% script loads image and runs a filter of a disc on them using 4 different
%%methods: imfilter conv and imfilter corr, xcorr2, filter2 (which uses
%%conv2)
%% like xcorr2 because it is larger (though probably can make the other
%% same size)
%% load data
F = 2;% total number of frames
% R = 457; C = 535;
R = 172; C = 130; %% flip R and C from TillVision
fname =  'G:\Scanziani Lab data (not Backup)\image data\022407_014\bksub\df2_47.tif' ;% name of first file
% df =
im = double(loadTILLimg(fname,C,R,F));
figure(1)
clf
imagesc(im);
colormap('gray')
%% make disc
D = makedisc(5);
% D = uint16(~D).*ones(size(D),'uint16').*uint16(max(max(im))); % invert
% D = D.*2^16;%uint16(max(max(im))); % invert
% D(find(D)==0) = min(min(im));
figure(2)
clf
imagesc(D)
colormap('gray')
%%
temp =length(D);
% tic
% F1 = imfilter(double(im),double(D));
% temp1 = 'imfilter'
% toc
% figure(3)
% clf
% F1 = F1((temp):(end-temp),(temp):(end-temp));
% imagesc(F1);
% colormap('gray')
% 
% tic
% F2 = imfilter(double(im),double(D),'conv');
% temp1 = 'imfiltercv'
% toc
% figure(4)
% F2 = F2((temp):(end-temp),(temp):(end-temp));
% imagesc(F2)
% colormap('gray')
% tic

F3 = xcorr2(double(im),double(D));
temp1 = 'xcorr2';
% toc
figure(5)
% imagesc(F3(round(length(D)/2):(end-round(length(D)/2)),round(length(D)/2):(end-round(length(D)/2))))
F3 = F3((temp):(end-temp),(temp):(end-temp));
whos F3
imagesc(F3)
colormap('gray')

% tic
% F4 = filter2(double(D),double(im));
% temp1 = 'filter2'
% toc
% 
% figure(6)
% F4 = F4((temp):(end-temp),(temp):(end-temp));
% imagesc(F4)
% % caxis([(max(caxis)-(max(caxis)-min(caxis))/5) max(caxis)])
% colormap('gray')
% figure(7)
% hist(F4(1:end),256)
% 
% size(im)-length(D)
%%
% Threshold image.
imth = mean(mean(F3))+ mean(std(F3))*3
bF3 = F3>imth;
figure(5);imagesc(bF3);
% NOTE: threshold.m should work as a gui to impliment this but needs some
% modification so that the image can be directly passed to it.
%% Get objects
[B,L,N] = bwboundaries(bF3,'noholes');

% remove polygons that have perimeter < minPer)
minPer = 10; sel = [];i=1;
for k = 1:length(B)
    if length(B{k})> minPer
        selB{i} = B{k};
        i=i+1;
    end

end

% display polygons for user to check validity of selected cells    
colors=['b' 'g' 'r' 'c' 'm' 'y'];
for k=1:length(selB)
    boundary = selB{k};
        cidx = mod(k,length(colors))+1;
        hold on;
    plot(boundary(:,2), boundary(:,1),...
        colors(cidx), 'LineWidth',2);
hold all
    rndRow = ceil(length(boundary)/(mod(rand*k,7)+1));
    col = boundary(rndRow,2); row = boundary(rndRow,1);
    h = text(col+1, row-1, num2str(k));
    set(h,'Color',colors(cidx),...
        'FontSize',14,'FontWeight','bold');
end

%save
cellBoundary = selB;
savefile = 'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Data Analysis\Expts\Calcium Imaging\022407_014\092007 getrois\cellBnd.m';
save(savefile,cellBoundary);

% convert polygons to mask
for i = 1:length(selB)
    cell_mask(i,:,:) = poly2mask(selB{i}(:,2),selB{i}(:,1),size(F3,1),size(F3,2));
end
