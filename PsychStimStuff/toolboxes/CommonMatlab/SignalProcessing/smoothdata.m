function [df dsmooth] = smoothdata(df,smwind)
% function [df dsmooth] = smoothdata(df,smwind)
% smooth data with multiple sweeps 
% INPUT: df matrix with sweep in each row

if min(size(df))==1 & ~isrowvector(df)
    df = df';
end
dsmooth = zeros(size(df));
for i =1:size(df,1)
    dsmooth(i,:) = smooth(df(i,:),smwind,'lowess')';
end
df = df -dsmooth ;

