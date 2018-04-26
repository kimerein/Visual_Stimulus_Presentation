%%
% plotsummaryofEPSCIPSCs.m

%%Get index of peaks that 
w= 4;
%% make event time monotonicsly increasing rather then restart at each
%% sweep
lengthsw = 4000 ;%% USER ENTRY %% this should be the length of a sweep in ms
temp1 = event(A).time+(event(A).Trace-1).*lengthsw;
temp2 = event(B).time+(event(B).Trace-1).*lengthsw;
ind = findSimEvent(temp1,temp2,w);
N= length(ind);
%% 
%% Subset of data
%% find half maximum amplitude 
%% to compare whethere there is a gerater correlation at high or low
%% amplitudes
% xtemp = A; temp = max(abs(event(xtemp).peak(ind(:,1))))/2*    max(event(xtemp).peak(ind(:,1)))/abs(max(event(xtemp).peak(ind(:,1))));
% temp = -80
% sel_ind = find(event(A).peak(ind(:,1))<temp);
% x = event(A).peak(ind(sel_ind,1));
% y = event(B).peak(ind(sel_ind,2));

%% ALL DATA
x = event(A).peak(ind(:,1));
y = event(B).peak(ind(:,2));
p = polyfit(x,y,1);
xf = [0:1/10*max(x)/abs(max(x)):max(abs(x))*max(x)/abs(max(x))]';
f = polyval(p,xf);;
figure(1)
clf
plot(xf,f,'-b')
hold on;
plot(x,y,'.b');
%% unity
% line([0 600],[0 600],'Color','k')
xlabel('IN0 Current Amp (pA)')
ylabel('IN1 Current Amp (pA)')
% ylabel('IN1 Ex Current Amp (pA)')
% xlabel('IN0 Ex Current Amp (pA)')
conf = 0.01;
r = corrcoef(x,y,'alpha',conf);
title(['y = ' num2str(p(1),'%1.2f') 'x ' num2str(p(2),'%1.2f') ' R:' num2str(r(2,1),'%1.2f') ' (p<' num2str(conf,'%1.2f') ') N=' num2str(N,'%d')],'Interpreter','none');

savetype = 'emf';
figString = ['SynScal_' event(A).expnum '_' event(A).sfilename '_' event(B).sfilename ]
if breport
    dirtemp = 'Pairs';
    figdesc = [figString];
    savefigure(writedirheader,dirtemp,figdesc,savetype)
end
