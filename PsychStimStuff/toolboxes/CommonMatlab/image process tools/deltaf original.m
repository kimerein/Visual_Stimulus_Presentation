function seq_out = deltaf(seq_in,startavg,numberavg,scale);
% function seq_out = deltaf(seq_in,startavg,numberavg,scale);
% Calculates  ((F - AvgF)/AvgF +1)*scale
% INPUT
%       seq_in: RxCxF, 
%       startavt: first frames F to average
%       numberavg: number of frames F to average
%       scale: factor to scale by
% 
% 2007_02_08 Sebi
% BA020907 modified

seq_out = zeros(size(seq_in,1),size(seq_in,2),size(seq_in,3),'int16'); % initialize output matrix
f_mean  = zeros(size(seq_in,1),size(seq_in,2),'int16');                % initialize mean matrix

% sum up 'numberavg' images starting with frame 'startavg'
for i = 1:numberavg
    f_mean = f_mean + seq_in(:,:,startavg+i-1);
end

% divide sum by the number of images 'numberavg'
f_mean = (f_mean/numberavg);
% f_mean = mean(seq_in(:,:,startavg:numberavg-1,3);

% calculate output pixelwise for every frame and sclae output
for i = 1: size(seq_in,3)  %% slower but less memory then matrix operation
    seq_out(:,:,i) = scale * ( (seq_in(:,:,i)-f_mean)./f_mean + 1);
end
    