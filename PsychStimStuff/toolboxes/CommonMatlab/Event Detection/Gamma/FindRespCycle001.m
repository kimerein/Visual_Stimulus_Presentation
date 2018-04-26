load 'a6831pir_in_vivo 006 -80mV odor1_MAT' Resp Current;
Ts = 5e3;
Nsamples = length(Resp);

%%% HIGH PASS FILTER
HIGHCUT = 1;%Hz
[B,A] = butter(2,2*HIGHCUT/Ts,'high'); 
% Y = filtfilt(B,A,Resp(1,:));
Y = filtfilt(B,A,Resp);
Y = Y';
Current = Current';
%%% LOW PASS FILTER
LOWCUT = 10;%Hz
[B,A] = butter(2,2*LOWCUT/Ts,'low'); 
Y1 = filtfilt(B,A,Y);  
%% Y1 is used for detecting breathing cycle by threshold only

LOWCUT = 300;%Hz
[B,A] = butter(2,2*LOWCUT/Ts,'low'); 
Y = filtfilt(B,A,Y);
% Y is used for aligning breathing cycles

% figure(13)
% plot(Y)

% figure(10);
% clf
% plot(Y-5*max(Y))
% hold on 
% plot(Resp(1,:),'r');

rthres = .001; %% USER DEFINE per experiment
j=1;
index = find(diff((Y1>rthres))==1)  %% find the index in the array of when signal crosses threshold

clear Y1;

bdebug = 0;

%% debug
if bdebug
    figure(11)
    clf
    plot(index,rthres.*ones(1,length(index)),'.r')
    hold on
    plot(Y)
    axis tight
end


%% FIND PEAK for each threshold crossing
ResWOI = [-.2 .5].*Ts; %% window of interest in index
PkWOI = [0 .1].*Ts; %%window to look for peak
                   %% code assumes PkWOI is small then ResWOI
clear Breath; clear indPeak; ii = 0; skipped = [];
for i=1:length(index)
    %%% ADD FIND PEAK 
        %% this point is used because it seems reproducible and high SNR in
        %% this data set.  A better measure might be slope 
        
    %% CHECK ResWOI exists in the same sweep
    %% index is a linear index from Resp which has a new Row for each sWeep
    %% so must find the linear index that corresponds to the edges of
    %% current sweep
    Nsweep = floor(index(i)/Nsamples) + 1;
    %% PROB .. loosing high end peaks
    
    if((index(i)+ResWOI(1))>(Nsweep-1)*Nsamples) && (index(i)+ResWOI(2))<Nsamples*Nsweep
                ii=ii+1;
        indPeak(ii) = find(Y(index(i)-PkWOI(1):(index(i)+PkWOI(2)))==max(Y(index(i)-PkWOI(1):(index(i)+PkWOI(2)))),1,'first') + index(i)-PkWOI(1)-1
            %% Get index of peak

        Breath(ii,:) = Y((indPeak(ii)+ResWOI(1)):(indPeak(ii)+ResWOI(2))); %% collect waveform of breath
        BCurr(ii,:) = Current((indPeak(ii)+ResWOI(1)):(indPeak(ii)+ResWOI(2)));
    else
        skipped = [skipped ; index(i) -1]; % code for skipped because of no ResWOI
    end
end

%% debug
if bdebug
    figure(11)
    clf
    plot(indPeak,rthres.*ones(1,length(indPeak)),'.r')
    hold on
    plot(skipped(:,1),rthres.*ones(1,size(skipped,2)),'.g')
    hold on
    plot(Y)
    axis tight
    figure(12)
    plot(Breath')
    axis tight
    
end


%%%%%%%%%%%%%%%%%%%%%
% APPLY criteria to remove weird breathing cycle


%%METHOD: Find 
indmax = find(RT == max(RT));

intBreath = sum(abs(Breath')) %% collapse breathing waveform into integral of absolute value
indexWeird = find(intBreath > 2*min(intBreath))  %% Note assumes that weird breathing is always larger in amplitude
                                                %% and regular breathing is
                                                %% all very similar to
                                                %% minium

if ~isempty(indexWeird)
    rejected = [indexWeird' -2*ones(length(indexWeird),1)]
end

temp = [1:size(Breath,1)];
temp(indexWeird)=0;
accepted = find(temp>0);


    
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




