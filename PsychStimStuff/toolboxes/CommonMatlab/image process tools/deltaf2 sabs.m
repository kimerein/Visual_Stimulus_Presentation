function seq_out = deltaf2(seq_in,startavg,numberavg,scale);

% no background subtraction included 
% 
% 2007_02_08 Sebi
seq_out = zeros(size(seq_in,1),size(seq_in,2),size(seq_in,3)); % initialize output matrix
f_mean  = zeros(size(seq_in,1),size(seq_in,2));                % initialize mean matrix

% sum up 'numberavg' images starting with frame 'startavg'
for i = 1:numberavg
    f_mean = f_mean + seq_in(:,:,startavg+i-1);
end

% divide sum by the number of images 'numberavg'
f_mean = (f_mean/numberavg);

% calculate output pixelwise for every frame and sclae output
for i = 1: size(seq_in,3)  
    seq_out(:,:,i) = scale * ( (seq_in(:,:,i)-f_mean)./f_mean + 1);
end
    