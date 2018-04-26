%% plot LFP and LFP triggered Ex and In currents
%% plot lag of Ex and In
%% requires summarydata to be loaded from .mat file
figure(301)
clf
ii=1; jj=1;
clear exlag; clear inlag
for i = 1:size(summarydata,2)-1
%  i=16
summarydata(i).exptnum
    if summarydata(1,i).holdVoltage > -50
        temptmean = -summarydata(i).tmean;
        xtime = [1:length(temptmean)].*summarydata(i).dt*1000;
        temptmean = temptmean-mean(temptmean);
        emean = summarydata(i).emean;
        emean = emean - mean(emean);
        emean = emean/max(emean);
        temptmean  = temptmean /max(temptmean)
        plot(xtime,temptmean,'-b','LineWidth',2)
        plot(xtime,temptmean,'-b','LineWidth',2)
        hold on
        inlag(ii,:) = [summarydata(1,i).tmeanlag(1)  summarydata(1,i).tmeanlagpk(1)];
        ii = ii+1
    else
        temptmean = summarydata(i).tmean;

        temptmean = temptmean-mean(temptmean);
        temptmean  = temptmean /max(temptmean);
        xtime = [1:length(temptmean)].*summarydata(i).dt*1000;

        plot(xtime,temptmean,'-r','LineWidth',2)
        hold on

        exlag(jj,:) = [summarydata(1,i).tmeanlag(1)  summarydata(1,i).tmeanlagpk(1)];
        jj = jj +1;
        pause
    end
end
mIN = [mean(inlag) sem(inlag)]
mEX = [mean(exlag) sem(exlag)]

figure;
plot([1 1], [mIN(1,1) mEX(1,1)],'.b','MarkerSize',20);
hold on
errorbar(1, mIN(1,1),mIN(1,3))

errorbar(1, mEX(1,1),mEX(1,3))
ylim([0 10])
xlim([.5 1.5])

figure;
plot([1 1], [mIN(1,2) mEX(1,2)],'.b','MarkerSize',20);
hold on
errorbar(1, mIN(1,2),mIN(1,4))

errorbar(1, mEX(1,2),mEX(1,4))
ylim([0 10])
xlim([.5 1.5])


clear EXisi;
clear INisi;
i=1;j=1;
for z = 1:8
    % clear msumdata;
    ind = find(~isnan(event(z).isi));
    x = event(z).isi(ind);
    y = event(z).peak(ind);
    y = circshift(y,1);
    N = length(x);
    p = polyfit(x,y,1);
    xf = [0:1/10*max(x)/abs(max(x)):max(abs(x))*max(x)/abs(max(x))]';
    xf = [0:1/10*max(x)/abs(max(x)):max(abs(x))*max(x)/abs(max(x))]';
    f = polyval(p,xf);;
    figure(3)
    clf
    plot(xf,f,'-r')
    hold on;
    plot(x,y,'.r');
    %% unity
    % line([0 600],[0 600],'Color','k')
    xlabel('ioi (ms)')
    ylabel('pk E/IPSC(pA)')
    % ylabel('IN1 Ex Current Amp (pA)')
    % xlabel('IN0 Ex Current Amp (pA)')
    conf = 0.01;
    r = corrcoef(x,y,'alpha',conf);
    title(['y = ' num2str(p(1),'%1.2f') 'x ' num2str(p(2),'%1.2f') ' R:' num2str(r(2,1),'%1.2f') ' (p<' num2str(conf,'%1.2f') ') N=' num2str(N,'%d')],'Interpreter','none');

    if event(z).V <-50
        EXisi(i,:) = [r(2,1) p(1)];
        i=size(EXisi,1)+1;

    else
        INisi(j,:) = [r(2,1) p(1)];
        j=size(INisi,1)+1;
    end
end
% 
% A = 2;
% B = 6;
% event(A).V
% event(B).V
% msumdata(i,:) =[mean(event(A).peak) mean(event(B).peak) std(event(A).peak) std(event(B).peak)];
% i=size(msumdata,1)+1;