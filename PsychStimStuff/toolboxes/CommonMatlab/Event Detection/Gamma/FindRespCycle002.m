clear all;
load 'a6831pir_in_vivo 006 -80mV odor1_MAT' Resp Current;
Ts = 5e3;
Nsamples = length(Resp);







%%% PRE PROCESSING Breathing for thresholding

%%% HIGH PASS FILTER
HIGHCUT = 1;%Hz
[B,A] = butter(2,2*HIGHCUT/Ts,'high');
%%% LOW PASS FILTER
LOWCUT = 10;%Hz for Y1
[BL,AL] = butter(2,2*LOWCUT/Ts,'low');
LOWCUT = 300;%Hz  for Y
[BLL,ALL] = butter(2,2*LOWCUT/Ts,'low');
for i=1:size(Resp,1)
    Y(i,:) = filtfilt(B,A,Resp(i,:));
    Y1(i,:)= filtfilt(BL,AL, Y(i,:));%% Y1 is used for detecting breathing cycle by threshold only
    Y(i,:)= filtfilt(BLL,ALL, Y(i,:));% Y is used for aligning breathing cycles
end
% Y = filtfilt(B,A,Resp(1,:));

Y = Y';
Y1 = Y1';
Current = Current';


%%%%%%%% THRESHOLD BREATHING CYCLES

rthres = .001; %% USER DEFINE per experiment
j=1;
index = find(diff((Y1>rthres))==1) ; %% find the index in the array of when signal crosses threshold

clear Y1;

bdebug = 1;


%% FIND PEAK for each threshold crossing
ResWOI = [-.05 .5].*Ts; %% window of interest in index
PkWOI = [0 .1].*Ts; %%window to look for peak
%% code assumes PkWOI is small then ResWOI
clear Breath; clear indPeak; ii = 0; skipped = [];
tic
for i=1:length(index)

    Nsweep = floor(index(i)/Nsamples) + 1;

    if((index(i)+ResWOI(1))>(Nsweep-1)*Nsamples) && (index(i)+ResWOI(2))<Nsamples*Nsweep
        %% CHECK ResWOI exists in the same sweep
        %% index is a linear index from Resp which has a new Row for each sWeep
        %% so must find the linear index that corresponds to the edges of
        %% current sweep
        ii=ii+1;
        indPeak(ii) = find(Y(index(i)-PkWOI(1):(index(i)+PkWOI(2)))==max(Y(index(i)-PkWOI(1):(index(i)+PkWOI(2)))),1,'first') + index(i)-PkWOI(1)-1;
        %%% FIND PEAK
        %% this point is used because it seems reproducible and high SNR in
        %% this data set.  A better measure might be slope %% Get index of peak

        Breath(ii,:) = Y((indPeak(ii)+ResWOI(1)):(indPeak(ii)+ResWOI(2))); %% collect waveform of breath
        BCurr(ii,:) = Current((indPeak(ii)+ResWOI(1)):(indPeak(ii)+ResWOI(2)));
        TimeinSweep(ii) = indPeak(ii) - Nsamples*(Nsweep-1);
    else
        skipped = [skipped ; index(i) -1]; % code for skipped because of no ResWOI
    end
    
%     if mod(i,20)==0
%         i
%         toc
%         tic
%     end
end

%% debug
if bdebug
    figure(11)
    clf
    plot(indPeak,rthres.*ones(1,length(indPeak)),'.r')
    title(['NumB: ' num2str(length(indPeak)) ' NumBsk: ' num2str(size(skipped,2))]);
    hold on
    plot(skipped(:,1),rthres.*ones(1,size(skipped,1)),'.g')
    hold on
    plot(Y(1:end))
    axis tight
    figure(12)
    plot(Breath')
    axis tight

end


%%%%%%%%%%%%%%%%%%%%%
% REMOVE weird breathing cycle
intBreath = sum(abs(Breath')); %% collapse breathing waveform into integral of absolute value
indexWeird = find(intBreath > 2*min(intBreath));  %% Note assumes that weird breathing is always larger in amplitude
%% and regular breathing is
%% all very similar to
%% minium

if ~isempty(indexWeird)
    rejected = [indexWeird' -2*ones(length(indexWeird),1)];
end

temp = [1:size(Breath,1)];
temp(indexWeird)=0;
accepted = find(temp>0); %% index in Breath of accepted Breaths


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%USER DEFINE  when ODOR comes on
odorStimTime = Ts*[1 3];
%% describe breathCat columns

% for i=1:size(odorTime,1)
breathCat = zeros(length(indPeak),size(odorStimTime,1)*size(odorStimTime,2),'int16');

breathCat(:,1)= TimeinSweep<odorStimTime(1,1); %% baseline
breathCat(:,2)=  TimeinSweep>=odorStimTime(1,1)& TimeinSweep<odorStimTime(1,2); %% during odor
breathCat(:,3)=  TimeinSweep>=odorStimTime(1,2)& TimeinSweep<(odorStimTime(1,2)+2*Ts); %after odor

%%removed rejected breaths
breathCat(indexWeird,:) = 0;

%% CHeck breaths are correctly cata
ind_tmp = find(breathCat(:,1)); %%
figure(102)  %% plot breaths that are catagorized as ind_tmp
for(i=1:size(Y,2))
    %%% FIX!!!! 
    plot(TimeinSweep(ind_tmp),rthres.*ones(1,length(TimeinSweep(ind_tmp)))-(i-1)*5e-3.*ones(1,length(TimeinSweep(ind_tmp))),'.r')
    hold on
    plot(Y(:,i)-(i-1)*5e-3)
    axis tight
end

%% 
figure(101)
clf
ind_tmp = find(breathCat(:,1)); %%  GET THE INDEXES  to plot
plot([1:length(Breath)]/Ts,BCurr(ind_tmp,:),'-k');
hold on
mBCurrBase = mean(BCurr(ind_tmp,:),1);
sBCurrBase = std(BCurr(ind_tmp,:),1);



ind_tmp = find(breathCat(:,2)); %%
plot([1:length(Breath)]/Ts,BCurr(ind_tmp,:),'-k');
mBCurr = mean(BCurr(ind_tmp,:),1);
sBCurr = std(BCurr(ind_tmp,:),1);
plot([1:length(BCurr)]/Ts,mBCurr,'-r','LineWidth',2);

plot([1:length(BCurr)]/Ts,mBCurrBase,'-b','LineWidth',2);
axis tight
mBreath = mean(Breath(ind_tmp,:),1);
plot([1:length(BCurr)]/Ts,mBreath/max(mBreath)*max(mBCurr)+max(mBCurr),'-g','LineWidth',1);


% plot([1:length(BCurr)]/Ts,mBCurr-sBCurr,'--k','LineWidth',1);




%%%%%%%%%%% plot breaths
figure (100)
clf
set(gcf,'Position',[700 0 1500 1120]);
subplot(2,1,1)
plot([1:length(Breath)]/Ts,Breath(rejected(:,1),:),'-k');
hold on
plot([1:length(Breath)]/Ts,Breath(accepted,:),'-b');
axis tight
subplot(2,1,2)
plot([1:length(BCurr)]/Ts,BCurr(accepted,:),'-b');
hold on
mBCurr = mean(BCurr(accepted,:),1);
sBCurr = std(BCurr(accepted,:),1);
plot([1:length(BCurr)]/Ts,mBCurr,'-r','LineWidth',2);
plot([1:length(BCurr)]/Ts,mBCurr-sBCurr,'--k','LineWidth',1);
axis tight




