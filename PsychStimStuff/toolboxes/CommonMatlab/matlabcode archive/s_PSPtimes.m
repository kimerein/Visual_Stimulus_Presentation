%% use with s_SpikeSort... 
%% after s_plotklusterspikes.m (where temp_spikeData is defined)

%%% Extract time of IPSC/ EPSC after spike event.
%% spikeData = time of spike
%% intraData = intracellular signal that PSP should be extracted from
%% 
bdebug =1;
bIPSC =0;
window = 5e-3; % window in Sec after spike to use for finding PSP max
X = 0.5 % 50% of max/min after spike

temp1 = zeros(size(temp_spikeData{j},1),int32(Ts*window)+1);
temp = output.data(:,2:output.Nchan:end);
% temp = temp - int16(mean(mean(temp,1)));%% not great way to get rid of the mean if there is variability in holding current or rest voltage
for jjjj=1:size(temp_spikeData{j},1)
    temp2 = (temp_spikeData{j}(jjjj,1));
    if temp2
        temp1(jjjj,:) = temp(temp2:(temp2 + int32(Ts*window)));
        temp1(jjjj,:) = smooth(temp1(jjjj,:),round(size(temp1,2)/13),'rloess');
    end

    %% oversample to get accurate values for parameters
    if jjjj ==1
        x = [1:size(temp1,2)];
        step = .1; xs = [1:step:size(temp1,2)];
    end
    data = pchip(x,temp1(jjjj,:),xs); %% if y is array size(y,2) == size(x)
    %% find time of 10% of max
    if bIPSC
        data = data -min(data); %% insure data >=0
        data = abs(data);

        PSP_time(jjjj) = find(abs(data -X*max(data)) == min(abs(data - X*max(data))));
    else
        data = data -max(data);

        PSP_time(jjjj) = find(abs(data - X*min(data)) == min(abs(data - X*min(data))));
    end

    if bdebug
        plot(data);
        line([PSP_time(jjjj) PSP_time(jjjj)],[min(data)  max(data)])
        %     line([1000 PSP_time(jjjj)],[min(data)  min(data)])
        title(num2str( PSP_time(jjjj) ))
        pause;
    end
    PSP_time = PSP_time*step *output.dt*1000;  %% in ms
end;