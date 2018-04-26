function [S,t,f] = plotcohgramCHX(df,df2,dt,movingwin,sparams,blogplot)
% function [S,t,f] = plotcohgramCHX(df,df2,dt,movingwin,sparams,blogplot)
% plot cohgram using Chronux cohgramc.
% INPUT: df, (each row contains new trial/sweep) vector/matrix 
%            NOTE: if multiple sweeps exist they are concatinated into one
%            long spectragram
%        dt sampling interval (usually in seconds)
%        movingwin = [windowsize step size ] (same units as dt)
%        sparams (optional) see mtspecgramc documentation
%
% BA061308
sparams.Fs = 1/dt; 
s = 'n';
if nargin >5 
    if blogplot;        s = 'l'; end
end

if min(size(df))==1 & ~isrowvector(df)
    df = df';
end
if min(size(df2))==1 & ~isrowvector(df2)
    df2 = df2';
end

[S,temp,temp,temp,temp,t,f]=cohgramc(df(:,:)',df2(:,:)',movingwin,sparams);

%plot spectragram
 if length(size(S))==2; tempS = S;xtime = [1:prod(size(df))]*dt;
 else [tempS t] = reduceDim(S,t); xtime = [1:length(df)]*dt;%% just for plotting not for 
 end
plot_matrix(tempS,t,f,s)
title('Coherogram')
xlabel('time(s)')
ylabel('frequency (Hz)')