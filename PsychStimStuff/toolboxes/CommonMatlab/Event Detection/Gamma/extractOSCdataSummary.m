DefineDir
expttype = 'KinateOsc';
writedirpath = [DATAANAL_DIR expttype '\'];
sumdata = 'SummaryData.mat';
load([writedirpath sumdata]);


%% find all same cell
dataindex = [];
% celldate = '2s6c1_2006_07_28';
breport = 1;
%         celldate = '2s5c1_2006_07_28';
        celldate = 's1c1_2_2006_08_09';
%         date = '2006_07_8';
for i=1:size(summarydata,2)  %% search for duplicates of the same analysis
    if strfind(summarydata(i).exptnum,celldate)
        dataindex = [dataindex i];
    end
end

clear title;
clear sexpt; clear data; clear Stmean; clear Semean; clear Sv; clear Slags;
figure(99)
clf
%% extract selected experiment data
for i=1:length(dataindex)
    in = dataindex(i);
    Stmean(:,i) = summarydata(in).tmean;
    Semean(:,i) = summarydata(in).emean;
    Sv(i) = summarydata(in).holdVoltage;
    sexpt{i} = summarydata(in).exptnum;
    %% lags
    Slags(i,:) = [summarydata(in).holdVoltage summarydata(in).tmeanlag]
end


%% MAKE PLOTS
clear data2; clear autoL; clear autoI;
for i=1:length(dataindex)
    if i ==1
        figure(200);clf; %% rawdata
        figure(201);clf; %% summary
    end
    %% SELECT color for holding
    in = dataindex(i);
    if summarydata(in).VC
        if summarydata(in).holdVoltage < -70
            xcolor = '-r';
        elseif summarydata(in).holdVoltage > -20
            xcolor = '-g';
        else
            xcolor = '-c';
        end
    else
        xcolor = '-k';
    end

    %% LFP and current
    figure(200)
    set(gcf,'Position',[580 50 1000 700],'PaperUnits','normalized');    hold on;    orient tall
    subplot(length(dataindex),2,(i-1)*2 +1)
    plot([1:size(summarydata(in).prodata1sec,2)].*output.dt*1000,summarydata(in).prodata1sec,'LineWidth',2)
    ylim([min(summarydata(in).prodata1sec(.1*1/output.dt:1/output.dt)) max(summarydata(in).prodata1sec(.1*1/output.dt:1/output.dt))])
    xlim([100 1000]);
    curylim1(i,:) = ylim;
    xlabel('time (ms)')
    title(['V:' num2str(Sv(i),'%1.1f')],'Interpreter','none')

    subplot(length(dataindex),2,i*2)
    temp = summarydata(in).intradata1sec - mean(summarydata(in).intradata1sec)
    plot([1:size(temp,2)].*summarydata(in).dt*1000,temp,xcolor,'LineWidth',2)
    ylim([min(temp(.1*1/summarydata(in).dt:1/summarydata(in).dt)) max(temp(.1*1/summarydata(in).dt:1/summarydata(in).dt))])
    xlim([100 1000]);
    curylim(i,:) = ylim;
    xlabel('Time (ms)')
    ylabel('Current (pA)')
    title([sexpt{i}],'Interpreter','none')

    figure(201);
    set(gcf,'Position',[0 50 1000 700],'PaperUnits','normalized');    hold on;    orient tall
    %% PLOT %% psd
    subplot(2,2,1)
    plot(summarydata(in).psdLFP(:,1),summarydata(in).psdLFP(:,2),xcolor,'LineWidth',2)
    xlabel('Hz');  ylabel('arb^2/Hz'); hold on; axis tight;
    xlim([0 100]);
    title('Power Spectral Density')
    %% PLOT  %% unnormalized Currents
    subplot(2,2,2)
    plot([1:size(Stmean,1)]*summarydata(dataindex(i)).dt.*1000,Stmean(:,i)-mean(Stmean(:,i)),xcolor,'LineWidth',2);
    hold on
    axis tight
    xlabel('ms')
    ylabel('pA')
    title('\gamma triggered <Current>')

    %% PLOT XCORR
    subplot(2,2,4)
    data2(:,i) = summarydata(in).xcorrLFPIntra;
    %     autoL(:,i) = summarydata(in).autocorrLFP;
    %     autoI(:,i) = summarydata(in).autocorrIntra;
    plot(([1:size(data2,1)]-((size(data2,1)-1)/2 +1)).*summarydata(dataindex(i)).dt.*1000,data2(:,i),xcolor,'LineWidth',2);
    hold all
    axis tight; xlabel('ms');
    %       xlim([-10 10])
    xlim([-30 30])
    title('xcorr LFP & \gamma<Current>')

    % PLOT %% coherence
    subplot(2,2,3)
    plot(summarydata(in).cohLFPIntra(:,1),summarydata(in).cohLFPIntra(:,2),xcolor,'LineWidth',2)
    xlabel('Hz');  ylabel('arb'); hold on; axis tight;
    xlim([0 100]);
        title('Coherence')
    %% PLOT %% coherence
%     subplot(2,2,4)
%     plot(summarydata(in).phaseLFPIntra(:,1),summarydata(in).phaseLFPIntra(:,2),xcolor,'LineWidth',2)
%     xlabel('Hz');  ylabel('Deg'); hold on; axis tight;
%     xlim([0 100]);
%         title('Phase')
% 
end

%% set y axis all the same scale on for raw currents
for i=1:length(dataindex)
    figure(200);
    subplot(length(dataindex),2,i*2);
    ylim([min(curylim(:,1)) max(curylim(:,2))]);
    subplot(length(dataindex),2,(i-1)*2+1);
    ylim([min(curylim1(:,1)) max(curylim1(:,2))]);
        curylim(i,:) = ylim;
end

figure(200);
figString = ['Summary001_' celldate];
if breport
    dirtemp = '';
    figdesc = [figString];
    savefigure(writedirpath,dirtemp,figdesc,savetype)
end
figure(201);
figString = ['Summary002_' celldate];
if breport
    dirtemp = '';
    figdesc = [figString];
    savefigure(writedirpath,dirtemp,figdesc,savetype)
end

% 
% %%%%%%%%%%%%%%%%%
% %% PLOT conductance
% %% assume we are at the reversal potential, so all the current is Ex/In
% %% ADD CODE to CORRECTION factor if not actually at reversal
% ExRev = 7.5; InRev = -87
% in1(1) = dataindex(find(v<-65,1)); in1(2) = dataindex(find(v> -30,1))
% is(1) = summarydata(in1(1)).holdVoltage- ExRev;
% is(2) = summarydata(in1(2)).holdVoltage - InRev;
% figure(100)
% clf
% for i=1:length(dataindex)
%     ind =find(in1== dataindex(i));
%     if(ind==1)
%         temp = Stmean(:,i)/is(1)
%         plot([1:size(Stmean,1)]*summarydata(dataindex(i)).dt.*1000,(temp)-mean(temp),'r')
%         hold all
%     end
%     if(ind==2)
%         temp = Stmean(:,i)/is(2)
%         plot([1:size(Stmean,1)]*summarydata(dataindex(i)).dt.*1000,temp-mean(temp),'g')
%         hold all
%     end
% end
% ylabel('\Delta nS') %% because there is ongoing synaptic input and it is hard to say when it is zero.
% xlabel('ms')
% title(['V_Ex :' num2str(ExRev) ' V_In :' num2str(InRev)]);
% 
% axis tight
