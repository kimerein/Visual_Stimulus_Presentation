function [stats fid] = plotFitCor(x,y,fid,bplot,cortype,plottype,pparam)
% function [stats fid] = plotFitCor(x,y,fid,bplot,cortype,plottype,pparam)
% INPUT: x and y are vectors, fid is scalar, if bplot is 0 doesn't plot, 
%  cortype (Pearson(1) or Spearman(0)). Pearson default 
%  plottype 1 hist3c else scatter (default scatter) 
%   dpath: DEFINEPATH for where plotset.m is located// CP 80115
%
% pparam.bins = x and y bin for hist3c plot <[30 30]>
% BA 112906

if nargin<5 || isempty(cortype)
    cortype = 1;
    if nargin <6 || isempty(plottype)
        plottype = 0;
    end
end

if nargin <7; pparam.bins = [30 30]; end   

N = length(x);
xf = [0:1/10*max(x)/abs(max(x)):max(abs(x))*max(x)/abs(max(x))]';

if size(x,1)==1
[b bint] = regress(y',[ones(length(x),1) x']);
else
[b bint] = regress(y,[ones(length(x),1) x]);
end
% p = polyfit(x,y,1);
p = circshift(b,1);
f = polyval(p,xf);
if nargin <4 
    bplot =1;
end
if ~cortype
    if isrowvector(x);x = x';end
    if isrowvector(y);y = y';end

    [r P] = corr(x,y,'type','Spearman','rows','pairwise');
else
    conf = 0.05;
[r P] = corrcoefnan(x,y,'alpha',conf);
end
if bplot
    if nargin <3 
        fid = figure;
    elseif isempty(fid)
        fid = figure;
    else
        figure(fid)
    end

    if plottype;        hist3c(x,y,pparam.bins,fid);
    else    plot(x,y,'ok');end
    hold on

    plot(xf,f,'-k')
    %   [R,P,RLO,RUP] = corrcoef(x,y,'alpha',conf);
    if cortype
        title(['y = ' num2str(p(1),'%1.2f') 'x +' num2str(p(2),'%1.2f') '  ,' num2str(bint(2,1),'%1.2f')  '<m<' num2str(bint(2,2),'%1.2f') '  R_P:' num2str(r(2,1),'%1.2f') ' (p<' num2str(P(2,1),'%1.2g') ') N=' num2str(N,'%d')]);
    else
    title(['y = ' num2str(p(1),'%1.2f') 'x +' num2str(p(2),'%1.2f') '  ,' num2str(bint(2,1),'%1.2f')  '<m<' num2str(bint(2,2),'%1.2f') '  R_S:' num2str(r,'%1.2f') ' (p<' num2str(P,'%1.2g') ') N=' num2str(N,'%d')],'Interpreter','none');
    end
    plotset(1,[]);
end

if cortype
    stats.r = r(2,1);
    stats.rP = P(2,1);
else
    stats.r = r;
    stats.rP = P;
end
stats.p = p;
stats.pconf = bint;
stats.N = N;

%% could alternatively boootstrap confidence interval
% bs = bootstrp(N, 'regress', y', [ones(length(X),1) X']);