% how to threshold image

% R = 457; C = 535;
R = 172; C = 130; %% flip R and C from TillVision
fname =  'G:\Scanziani Lab data (not Backup)\image data\022407_014\bksub\df2_47.tif' ;% name of first file
% df =
im = double(loadTILLimg(fname,C,R,F));
D = makedisc(5);
temp =length(D);

F3 = xcorr2(double(im),double(D));
temp1 = 'xcorr2'
figure(5)
F3 = F3((temp):(end-temp),(temp):(end-temp));
whos F3
imagesc(F3)
colormap('gray')

figure(6)
hist(F3(1:end),256)

th = mean(mean(F3))+std(std(F3))*7.5;
tF = zeros(size(F3));
tF(find(F3>th)) =1;
figure(7)
imagesc(tF)
BW_filled = imfill(tF,'holes');
[B L]= bwboundaries(BW_filled);
hold on
     for k = 1:length(B)
           boundary = B{k};
           plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
     end
     %% impose min  diameter of boundry
       
     
dim = size(tF)
col = round(dim(2)/2)-90
row = min(find(tF(:,col)))
%% 
% BW1 = edge(tF,'sobel');
BW2 = edge(tF,'canny');
figure(8)
imagesc(BW1)
figure, imshow(BW2)
BW = im2bw(tF);
figure(9)
imagesc(BW)
dim = size(BW)
col = round(dim(2)/2)-90;
row = min(find(BW(:,col)))