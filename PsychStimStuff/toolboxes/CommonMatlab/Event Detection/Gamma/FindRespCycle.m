
Ts = 5e3;
Nsamples = length(Resp);

%%% HIGH PASS FILTER
HIGHCUT = 1;%Hz
[B,A] = butter(2,2*HIGHCUT/Ts,'high'); 
Y = filtfilt(B,A,Resp(1,:));

Y = Y';
% figure(10);
% clf
% plot(Y-5*max(Y))
% hold on 
% plot(Resp(1,:),'r');

rthres = .003; %% USER DEFINE per experiment
j=1;
index = find(diff((Y>rthres))==1)  %% find the index in the array of when signal crosses threshold

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
    Nsweep = round(index(i)/Nsamples) + 1;
    %% PROB .. loosing high end peaks
    
    if((index(i)+ResWOI(1))>(Nsweep-1)*Nsamples) && (index(i)+ResWOI(2))<Nsamples*Nsweep
                ii=ii+1;
        indPeak(ii) = find(Y(index(i)-PkWOI(1):(index(i)+PkWOI(2)))==max(Y(index(i)-PkWOI(1):(index(i)+PkWOI(2)))),1,'first') + index(i)-PkWOI(1)-1
            %% Get index of peak

        Breath(ii,:) = Y((indPeak(ii)+ResWOI(1)):(indPeak(ii)+ResWOI(2))); %% collect waveform of breath
    else
        skipped(i) = -1 % code for skipped because of no ResWOI
    end
end

%% debug
if bdebug
    figure(11)
    clf
    plot(indPeak,rthres.*ones(1,length(indPeak)),'.r')
    hold on
    plot(Y)
    axis tight
end

%%%%%%%%%%%%%%%%%%%%%
% APPLY criteria to remove wierd breathing cycle
RT = Resp(1,17440:20570)  %% DEFINE template
indmax = find(RT == max(RT));



