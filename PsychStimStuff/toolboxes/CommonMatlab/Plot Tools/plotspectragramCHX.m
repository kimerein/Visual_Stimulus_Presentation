function [S,t,f] = plotspectragramCHX(df,dt,movingwin,sparams,blogplot,toffset)
% function [S,t,f] = plotspectragramCHX(df,dt,movingwin,sparams,blogplot,toffset)
% plot spectragram using Chronux mtspecgramc.
% INPUT: df, (each row contains new trial/sweep) vector/matrix 
%            NOTE: if multiple sweeps exist they are concatinated into one
%            long spectragram
%        dt sampling interval (usually in seconds)
%        movingwin = [windowsize step size ] (same units as dt)
%        sparams (optional) see mtspecgramc documentation
%        blogplot
%        toffset - offset in time (same units as dt) for ease in plotting
%        onto spectragram where data already exists
% BA061308
if ~exist('toffset','var');  toffset = 0; end
sparams.Fs = 1/dt; 
s = 'n';
if nargin >4 
    if blogplot
        s = 'l'; 
    end
end

[S,t,f]=mtspecgramc(df(:,:)',movingwin,sparams);

%plot spectragram
 if length(size(S))==2; tempS = S;xtime = [1:prod(size(df))]*dt;
 else [tempS t] = reduceDim(S,t); xtime = [1:length(df)]*dt;%% just for plotting not for 
 end
plot_matrix(tempS,t+toffset,f,s)
xlabel('time(s)')
ylabel('frequency (Hz)')