function [S,t,f] = plotspectragramsegCHX(df,dt,win,sparams,blogplot)
% function [S,t,f] = plotspectragramsegCHX(df,dt,win,sparams,blogplot)
% similar to plotspectragramCHX but segmented
%
% plot spectragram using Chronux mtspecgramsegc.
% INPUT: df, (each row contains new trial/sweep) vector/matrix 
%            NOTE: if multiple sweeps exist they are concatinated into one
%            long spectragram
%        dt sampling interval (usually in seconds)
%        win = time segment size (same units as dt)
%        sparams (optional) see mtspecgramsegc documentation
%
% BA112508

sparams.Fs = 1/dt; 
s = 'n';
if nargin >4 
    if blogplot
        s = 'l'; 
    end
end

[S,f]=mtspectrumsegc(df(:,:)',win,sparams,0);
t = ([1:size(S,2)]-1/2)*win; % time of each segment defined as middle segment
%plot spectragram
 if length(size(S))==2; tempS = S;
%  else [tempS t] = reduceDim(S,t); xtime = [1:length(df)]*dt;%% just for plotting not for 
 end
plot_matrix(tempS',t,f,s)
xlabel('time(s)')
ylabel('frequency (Hz)')