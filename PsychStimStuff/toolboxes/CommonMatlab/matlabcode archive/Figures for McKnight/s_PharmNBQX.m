%% NBQX spontaneous activity pharm
%%% CODE for MANUAL entering data
% 'a' is pasted from 500ms spiketimes NBQX.xls

% spiketimes_slice1_102105
N = 10; %secsperbin
data = bin2(a,N);  % Histogram raw spiketime data
% plot spike data with lines were pharmacology was applied
plot(data)

% hold
% clear sweeptime;
% M = 5 %sec per sweep
% sweepnum = [133 133+24 133+24+24]; %Slice1
% sweeptime = sweepnum*M/N;
% sweeptime = [sweeptime; ones(1,size(sweeptime,2))]';
% for i=1:size(sweepnum,2)
% line([sweeptime(i) sweeptime(i)],[1+1 max(data)],'Color','r','LineWidth',0.5); %VERTICAL
% end

% save('slice1_102105_003_004', 'spiketimes_slice1_102105','-mat');
% pharmdata = struct('Experiment','10_21_05_E1','DAMGO_baseline',[93*9-119 93*9],'DAMGO',[838 958],...
%     'CPA_baseline',[179*9-119 179*9],'CPA',[179*9+1 179*9+120], 'binned_spiketimes',spiketimes_slice1_102105);

%%% SET WINDOWS of baseline etc
% i = 2;
% pharmdata(i).Experiment = '022306';
% pharmdata(i).NBQX = [(36+24)*5 (36+48)*5];
% pharmdata(i).NBQX_baseline = [(36-24)*5 36*5];
% pharmdata(i).CPP_baseline = [(81-24)*5 81*5];
% pharmdata(i).CPP = [(81+24)*5 (81+48)*5];
% pharmdata(i).rawdata = single(a);
% pharmdata(i).binned_spiketimes = single(data); % 10sec per bin
% pharmdata(i).msperbin = N;

%% END OF DATA MANUAL ENTRY

%%  Normalize each experiment to the mean of the baseline
%% Print Windows for each data set, and
spikeDBase = zeros(size(pharmdata,2),(pharmdata(i).NBQX_baseline(1,2)-pharmdata(i).NBQX_baseline(1,1))/pharmdata(i).msperbin+1);
spikeD = zeros(size(pharmdata,2),(pharmdata(i).NBQX_baseline(1,2)-pharmdata(i).NBQX_baseline(1,1))/pharmdata(i).msperbin+1); %% 
selected = [1:2 4:9];
for i=1:size(pharmdata,2)
    if find(selected ==i)
    i
        temp1 = pharmdata(i).NBQX_baseline(1)/pharmdata(i).msperbin;
        temp2 = pharmdata(i).NBQX_baseline(2)/pharmdata(i).msperbin;
        tempmean = mean(pharmdata(i).binned_spiketimes(temp1:temp2));% normalize to mean activity
        spikeDBase(i,:) = pharmdata(i).binned_spiketimes(temp1:temp2)./tempmean;
        if isempty(pharmdata(i).CPP)
            temp1 = pharmdata(i).NBQX(1)/pharmdata(i).msperbin;
            temp2 = pharmdata(i).NBQX(2)/pharmdata(i).msperbin ;
        else
            temp1 = pharmdata(i).CPP(1)/pharmdata(i).msperbin;
            temp2 = pharmdata(i).CPP(2)/pharmdata(i).msperbin ;
        end
        spikeD(i,:) = pharmdata(i).binned_spiketimes((temp1:temp2))./tempmean;
%         figure(100)
%         plot([1:120],spikeDBase(i,1:120),[121:240],spikeD(i,1:120));
    
%     title(pharmdata(i).Experiment);
%         pause
    end
end

temp = spikeDBase;
temp1 = spikeD;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mDB = mean(temp);
sDB = sem(temp);
mD = mean(temp1);
sD = sem(temp1);

x = [1:size(mDB,2)];
x1 = [1:size(mD,2)]+max(x);
figure(201);
errorbar([x x1]*10,[mDB mD],[sDB sD],'.k',...
                'MarkerEdgeColor','b',...
                'MarkerFaceColor','b')
xlabel('time (s)')
axis tight

%%%mean Contrl and mean NBQX
mmDB = mean(mDB)
smDB = sem(mDB)
stmDB = std(mDB)
mmD = mean(mD)
smD = sem(mD)
stmD = std(mD)
figure(202);
bar([1 2],[mmDB mmD])
hold on
set(gca,'Color','none','XGrid', 'off','XTick',[1, 2],'XTickLabelMode','manual','XTickLabel',['Ctrl';'NBQX'])
errorbar([1 2],[mmDB mmD],[smDB smD],'.k',...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','b')
hold on;
axis tight


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% print out values
i =6;
pharmdata(i).Experiment 
pharmdata(i).DAMGO_baseline
pharmdata(i).DAMGO 
pharmdata(i).CPA_baseline
pharmdata(i).CPA 
pharmdata(i).DPCPX_baseline 
pharmdata(i).DPCPX 
pharmdata(i).NBQX_baseline =  [1116 1116+120]
pharmdata(i).NBQX = pharmdata(i).DAMGO_baseline
pharmdata(i).GABAZINE_baseline  = [204 324]
pharmdata(i).GABAZINE = [1116 1116+120]
pharmdata(i).TTX_baseline 
pharmdata(i).TTX 