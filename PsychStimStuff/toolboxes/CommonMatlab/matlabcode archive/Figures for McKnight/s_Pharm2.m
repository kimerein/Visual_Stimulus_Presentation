%% Slice Pharm for McKnight

% spiketimes_slice1_102105
N = 10; %secsperbin
data = bin2(a,N);  % Histogram raw spiketime data

% plot spike data with lines were pharmacology was applied
plot(data)
hold
clear sweeptime;
M = 5 %sec per sweep
sweepnum = [133 133+24 133+24+24]; %Slice1
sweeptime = sweepnum*M/N;
sweeptime = [sweeptime; ones(1,size(sweeptime,2))]';
for i=1:size(sweepnum,2)
line([sweeptime(i) sweeptime(i)],[1+1 max(data)],'Color','r','LineWidth',0.5); %VERTICAL
end

% save('slice1_102105_003_004', 'spiketimes_slice1_102105','-mat');
% pharmdata = struct('Experiment','10_21_05_E1','DAMGO_baseline',[93*9-119 93*9],'DAMGO',[838 958],...
%     'CPA_baseline',[179*9-119 179*9],'CPA',[179*9+1 179*9+120], 'binned_spiketimes',spiketimes_slice1_102105);

%%% SET WINDOWS of baseline etc
% i = 2;
% pharmdata(i).Experiment = '12_26_05_E1';
% pharmdata(i).DAMGO_baseline = [737 857];
% pharmdata(i).DAMGO = [858 978];
% pharmdata(i).CPA_baseline = [1060 1180];
% pharmdata(i).CPA = [1580-120 1580];
% pharmdata(i).DPCPX_baseline = [1704-120 1704];
% pharmdata(i).DPCPX = [2089-120 2089];
% pharmdata(i).NBQX_baseline = [580 700];
% pharmdata(i).NBQX = [1060 1180];
% pharmdata(i).GABAZINE_baseline = [10 130];
% pharmdata(i).GABAZINE = [180 300];
% pharmdata(i).TTX_baseline = [2641-120 2641];
% pharmdata(i).TTX = [2721 2843+57];
% pharmdata(i).rawdata = single(a);
% pharmdata(i).binned_spiketimes = single(data);


%% Normalize baseline (mean of baseline is 1)
%% Plot 
pharmdata(i).binned_spiketimes = single(pharmdata(i).binned_spiketimes);


%% Show windows chosen for baseline and data
%% (NOT necassary for the rest of the code .. just a display)
% for i=1:size(pharmdata,2)
%     i
%     %     Db = pharmdata(i).DAMGO_baseline(2) - pharmdata(i).DAMGO_baseline(1)
%     %     Db = pharmdata(i).DAMGO(2) - pharmdata(i).DAMGO(1)
%     %     Cb = pharmdata(i).CPA_baseline(2) - pharmdata(i).CPA_baseline(1)
%     %     Cb = pharmdata(i).CPA(2) - pharmdata(i).CPA(1)
%     figure(100)
%     plot(pharmdata(i).binned_spiketimes);
%     if ~isempty(pharmdata(i).NBQX_baseline)
%         line([pharmdata(i).NBQX_baseline(1) pharmdata(i).NBQX_baseline(1)],[1+1 max(pharmdata(i).binned_spiketimes)],'Color','r','LineWidth',0.5); %VERTICAL
%         line([pharmdata(i).NBQX_baseline(2) pharmdata(i).NBQX_baseline(2)],[1+1 max(pharmdata(i).binned_spiketimes)],'Color','r','LineWidth',0.5); %VERTICAL
%         line([pharmdata(i).NBQX(1) pharmdata(i).NBQX(1)],[1+1 max(pharmdata(i).binned_spiketimes)],'Color','c','LineWidth',0.5); %VERTICAL
%         line([pharmdata(i).NBQX(2) pharmdata(i).NBQX(2)],[1+1 max(pharmdata(i).binned_spiketimes)],'Color','c','LineWidth',0.5); %VERTICAL
%     end
%     if ~isempty(pharmdata(i).GABAZINE_baseline)
%         line([pharmdata(i).GABAZINE_baseline(1) pharmdata(i).GABAZINE_baseline(1)],[1+1 max(pharmdata(i).binned_spiketimes)],'Color','k','LineWidth',0.5); %VERTICAL
%         line([pharmdata(i).GABAZINE_baseline(2) pharmdata(i).GABAZINE_baseline(2)],[1+1 max(pharmdata(i).binned_spiketimes)],'Color','k','LineWidth',0.5); %VERTICAL
%         line([pharmdata(i).GABAZINE(1) pharmdata(i).GABAZINE(1)],[1+1 max(pharmdata(i).binned_spiketimes)],'Color','g','LineWidth',0.5); %VERTICAL
%         line([pharmdata(i).GABAZINE(2) pharmdata(i).GABAZINE(2)],[1+1 max(pharmdata(i).binned_spiketimes)],'Color','g','LineWidth',0.5); %VERTICAL
%         title(pharmdata(i).Experiment);
%     end
% 
%     pause;
% end

%%  Normalize each experiment to the mean of the baseline
%% Print Windows for each data set, and
spikeDBase = zeros(size(pharmdata,2),121);
spikeD = spikeDBase ;
spikeCPABase = spikeDBase ;
spikeCPA = spikeDBase ;
for i=1:size(pharmdata,2)
    i
    %     Db = pharmdata(i).DAMGO_baseline(2) - pharmdata(i).DAMGO_baseline(1)
    %     Db = pharmdata(i).DAMGO(2) - pharmdata(i).DAMGO(1)
    %     Cb = pharmdata(i).CPA_baseline(2) - pharmdata(i).CPA_baseline(1)
    %     Cb = pharmdata(i).CPA(2) - pharmdata(i).CPA(1)
    if ~isempty(pharmdata(i).DAMGO_baseline)
        tempmean = mean(pharmdata(i).binned_spiketimes(pharmdata(i).DAMGO_baseline(1):pharmdata(i).DAMGO_baseline(2)));% normalize to mean activity
        spikeDBase(i,:) = pharmdata(i).binned_spiketimes(pharmdata(i).DAMGO_baseline(1):pharmdata(i).DAMGO_baseline(2))./tempmean;
        spikeD(i,:) = pharmdata(i).binned_spiketimes(pharmdata(i).DAMGO(1):pharmdata(i).DAMGO(2))./tempmean;
%         figure(100)
%         plot([1:120],spikeDBase(i,1:120),[121:240],spikeD(i,1:120));
    end
    tempmean = mean(pharmdata(i).binned_spiketimes(pharmdata(i).CPA_baseline(1):pharmdata(i).CPA_baseline(2)));% normalize to mean activity
    spikeCPABase(i,:) = pharmdata(i).binned_spiketimes(pharmdata(i).CPA_baseline(1):pharmdata(i).CPA_baseline(2))./tempmean;
    spikeCPA(i,:) = pharmdata(i).binned_spiketimes(pharmdata(i).CPA(1):pharmdata(i).CPA(2))./tempmean;

%     figure(101)
%     plot([1:120],spikeCPABase(i,1:120),[121:240],spikeCPA(i,1:120));
%     title(pharmdata(i).Experiment);
%         pause
end

% Define experiment that should include DAMGO
temp = spikeDBase([1 3 5],:);
temp1 = spikeD([1 3 5],:);

% % MEAN changes in spikerate from baseline to drug
% mean(mean(spikeCPABase))
% std(std(spikeCPABase))
% mean(mean(spikeCPA))
% std(std(spikeCPA))
% meanDB = mean(mean(temp))
% std(std(temp))
% meanD = mean(mean(temp1))
% std(std(temp1))


%% Mean and SEM of all experiments and Bin into 10s bins
mCB = rebin(mean(spikeCPABase,1)',10);
sCB = rebin(sem(spikeCPABase,1)',10);
mC = rebin(mean(spikeCPA,1)',10);
sC = rebin(sem(spikeCPA,1)',10);
mDB = rebin(mean(temp)',10);
sDB = rebin(sem(temp)',10);
mD = rebin(mean(temp1)',10);
sD = rebin(sem(temp1)',10);

x = [1:size(mC,2)];
x1 = size(mC,2)+x;
figure(200);
errorbar([x x1],[mCB mC],[sCB sC],'.k',...
                'MarkerEdgeColor','b',...
                'MarkerFaceColor','b')
axis tight
figure(201);
errorbar([x x1],[mDB mD],[sDB sD],'.k',...
                'MarkerEdgeColor','b',...
                'MarkerFaceColor','b')
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