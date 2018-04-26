function PSP_time = PSPtime3(spiketime,intraData,dt,bIPSC,bIntra)
% function PSP_time = PSPtime(spiketime,intraData,dt,bIPSC,bIntra)
%% Extract time(ms) of IPSC/ EPSC after spike event.
%% spikeData = time of spike
%% intraData = intracellular signal that PSP should be extracted from
%% dt = time between samples
%% bIPSC = 1 if IPSCs
%% bIntra = set to 1 if intraData is not raw but already the PSP
%% fits the rise time with an exponential.
Ts = 1/dt;

HPFilter = 300;
[B,A] = butter(2,2*HPFilter/Ts,'low');
% intraData1 = intraData;
intraData = filtfilt(B,A,intraData);
HPFilter = 4;
[B,A] = butter(2,2*HPFilter/Ts,'high');
intraData = filtfilt(B,A,intraData);

PSP_time = [];
bdebug = 0;
 if nargin <4
    bIPSC =0;
end
window = 5e-3; % window in Sec after spike to use for finding PSP max
X = 0.5; % 50% of max/min after spike

if ~bIntra

    temp1 = zeros(size(spiketime,1),int32(Ts*window)+1);
else
    temp1 = zeros(size(spiketime,1),size(intraData,1));
end
% temp = temp - int16(mean(mean(temp,1)));%% not great way to get rid of the mean if there is variability in holding current or rest voltage
for jjjj=1:size(spiketime,1)
    temp2 = (spiketime(jjjj,1));
    if temp2
        if ~bIntra
            temp1(jjjj,:) = intraData(temp2 - int32(Ts*2*window):(temp2 + int32(Ts*window)));
        else
            temp1(jjjj,:) = intraData(:,jjjj);
        end
        %         temp1(jjjj,:) = smooth(temp1(jjjj,:),round(size(temp1,2)/13),'rloess');

    end
%     spike is at 10ms
% %% look 5ms before spike
%     ind_max = find(temp1(jjjj,int32(Ts*1*window):int32(Ts*2*window)+int32(Ts*2*window)+1) == max(temp1(jjjj,int32(Ts*1*window):int32(Ts*2*window)+int32(Ts*2*window)+1)),1,'first')+ int32(Ts*1*window) -1;
%     ind_min = find(temp1(jjjj,int32(Ts*1*window):int32(Ts*2*window)+int32(Ts*2*window)+1) == min(temp1(jjjj,int32(Ts*1*window):int32(Ts*2*window)+int32(Ts*2*window)+1)),1,'first')+ int32(Ts*1*window) -1;
   ind_max = find(temp1(jjjj,int32(Ts*2*window):int32(Ts*2*window)+int32(Ts*2*window)+1) == max(temp1(jjjj,int32(Ts*2*window):int32(Ts*2*window)+int32(Ts*2*window)+1)),1,'first')+ int32(Ts*2*window) -1;
    ind_min = find(temp1(jjjj,int32(Ts*2*window):int32(Ts*2*window)+int32(Ts*2*window)+1) == min(temp1(jjjj,int32(Ts*2*window):int32(Ts*2*window)+int32(Ts*2*window)+1)),1,'first')+ int32(Ts*2*window) -1;

    if     bIPSC
        PSP_time(jjjj) = ind_max;
    else
        PSP_time(jjjj) = ind_min;
    end
    %     PSP_time(jjjj) = ind_max;

    %     %% find time of 10% of max
    % %     if bIPSC
    %         temp1(jjjj,:) = temp1(jjjj,:) -min(temp1(jjjj,:))+1; %% insure data >=0
    %         temp1(jjjj,:) = abs(temp1(jjjj,:));
%         y = temp1(jjjj,1:ind_max)
%         x = [1:size(y,2)];
%         p = polyfit(x,log(y),1);
%         xf = [0:0.1:size(x,2)]';
%         f = polyval(p,xf);
%         plot(x,log(y),'ob',xf,f,'-g')
%     else
%         data = data -max(data);
%         ind_min = find(data == min(data));
%         PSP_time(jjjj) = find(abs(data(1:ind_min) - X*min(data)) == min(abs(data(1:ind_min) - X*min(data))));
%     end

    if bdebug
        figure(99);
        plot(temp1(jjjj,1:Ts*25e-3));
        line([PSP_time(jjjj) PSP_time(jjjj)],[min(temp1(jjjj,:))  max(temp1(jjjj,:))])
        %     line([1000 PSP_time(jjjj)],[min(data)  min(data)])
        title(num2str( PSP_time(jjjj) ))
%         pause;
    end
end;

% print('-dpdf', 'tempbass')
step =1;
    PSP_time = (PSP_time-(double(Ts*2*window) -1))*step *dt*1000;  %% in ms
