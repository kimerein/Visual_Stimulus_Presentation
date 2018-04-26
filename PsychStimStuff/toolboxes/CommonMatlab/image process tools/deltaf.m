function seq_out = deltaf(seq_in,startavg,numberavg,scale);
% function seq_out = deltaf(seq_in,startavg,numberavg,scale);
% Calculates  ((F - AvgF)/AvgF +1)*scale
% INPUT
%       seq_in: RxCxF, (uint16)
%       startavg: first frames F to average
%       numberavg: number of frames F to average
%       scale: factor to scale by
% 
% OUTPUT 
%       seq_out: uint16
% Note on data types and what scal factor to choose
%   images imported from tif (which is what this is made for) contain uint
%   format.  It is assumed that there are uint16. (there may be problems if
%   they are uint32 or geater.
%   Since some pixel values will be less then the mean, and dividing deltf
%   by f may result in fractions, the deltaF/F computation is done in
%   doubles (singles would probably suffice if it starts to be slow or too
%   large). All pixel values should be postive.  They won't necessarily be integers.  
%   Scale factor should be chosen so that the maximal possible resolution (which for 12-bit TILL is 1/4096?)
%   So an ideal scale factor would but something larger then 4096.
%   HOWEVER, for some reason when scaling data by great then 1000-5000 one
%   starts to get histograms (of intensity) that are are gaussian but have
%   spikes sticking out of them (maybe because of rounding problems between floating and integer).  Since for my experiments I do not get
%   close to the whole dynamic range (4096) I have settled on scale = 1000
%
% 2007_02_08 Sebi
% BA020907 modified

%% check uint16
%%

seq_out = zeros(size(seq_in,1),size(seq_in,2),size(seq_in,3),'uint16'); % initialize output matrix
f_mean  = zeros(size(seq_in,1),size(seq_in,2),'double');                % initialize mean matrix

% sum up 'numberavg' images starting with frame 'startavg'
for i = 1:numberavg
    f_mean = f_mean + double(seq_in(:,:,startavg+i-1));
end

% divide sum by the number of images 'numberavg'
f_mean = (f_mean/numberavg);
% f_mean = mean(seq_in(:,:,startavg:numberavg-1,3);

% calculate output pixelwise for every frame and sclae output
for i = 1: size(seq_in,3)  %% slower but less memory then matrix operation
    temp = scale * ( (double(seq_in(:,:,i))-f_mean)./f_mean + 1); % compute deltaf
%     %%check postive
%     if any(temp<1)
%         error('deltaf contains negative value (should be positive)');
%     end
     seq_out(:,:,i) = uint16(temp);
end
    