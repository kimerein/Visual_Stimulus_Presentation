function [Params Failure] = PSPParam(spiketime,intraData,dt,bIPSC,bIntra)
% function Params = PSPParam(spiketime,intraData,dt,bIPSC,bIntra)
%% Extract Parameters of PSP time(ms) of IPSC/ EPSC after spike event.
%% INPUT:
%%      spiketime is a vector of spiketimes in data set.
%%      intraData can have 2 forms
%%              if bIntra = 1
%%                  intraData is a matrix where EACH ROW contains
%%                  intracellular data from each spiketime.  spiketime
%%                  should occur 10ms into each row of intracellular data
%%              if bIntra = 0
%%                  intraData a matrix taken from intracellular data from abfstruct.data.
%%                  i.e. each column is a new episode/sweep 
%%       dt 
%%          sampling interval in secounds.
%%       bIPSC 
%%          set to 1 if EPSP or IPSCs are being detected otherwise set to 0
%%  OUTPUT:
%%       Params = 
%%             [ time of crossing of threshold
%%               time of peak] 
%%              cases no PSP:
%%                 NO baseline crossing are given params of -1/10
%%                 LOCAL peak between threshold-1ms and peak-0.25ms -2/10
%%
%%      Failure = spikes_ind where there was not voltage above/below
%%          threshold
bdebug = 0;
if nargin <4
    bIPSC =0;
end    
Ts = 1/dt;
Failure = [];
Params = zeros(size(spiketime,1),2,'single');
%% high pass filter
HPFilter = 4;
[B,A] = butter(2,2*HPFilter/Ts,'high');
intraData = filtfilt(B,A,intraData);

if bdebug
    intraData1 = intraData;
end
%% low pass filter
HPFilter = 400;
[B,A] = butter(2,2*HPFilter/Ts,'low');
intraData = filtfilt(B,A,intraData);


window = 10e-3; % window in Sec after spike to use for finding PSP max
 
if ~bIntra
    temp1 = zeros(size(spiketime,1),int32(Ts*window)+1);
else
    temp1 = zeros(size(spiketime,1),size(intraData,1));
end
% temp = temp - int16(mean(mean(temp,1)));%% not great way to get rid of the mean if there is variability in holding current or rest voltage
for jjjj=1:size(spiketime,1)
    temp2 = (spiketime(jjjj,1));
        % find baseline noise
    %% PARAMETERS
    baswindow =  10e-3; %define baseline window
    S = std(temp1(jjjj,1:int32(Ts*window)));
    Xstd = 3;
ind_thres =0;ind_peak =0;ind_HPK=0;
    if temp2
        if ~bIntra
            temp1(jjjj,:) = intraData(temp2 - int32(Ts*window):(temp2 + int32(Ts*window)));
            %% TODO catch if window is outside of what exists in data file
        else
            temp1(jjjj,:) = intraData(:,jjjj);
        end
    end
    temp1(jjjj,:) = temp1(jjjj,:)-mean(temp1(1:int32(Ts*window)));
  
% 
    %     ASSUMES spike is at 10ms into intraData
    try
        if     bIPSC %% NOTE: looks for threshold crossing from 0-10ms after spike
            % Detect PSP by Xstd from std.
            temp = find(temp1(jjjj, int32(Ts*window:Ts*2*window))> S*Xstd,1,'first');
            if (isempty(temp))
                Failure = [Failure jjjj];
                Params(jjjj,:) = -1/10;
            else
                ind_thres = int32(temp + Ts*window -1);
                peak = max(temp1(jjjj,ind_thres:int32(Ts*2*window)));
                ind_peak = find(temp1(jjjj,ind_thres:int32(Ts*2*window)) == peak,1,'first') + ind_thres-1;
            end
        else
            temp = find(temp1(jjjj, int32(Ts*window:Ts*2*window))< -S*Xstd,1,'first');
            if (isempty(temp))
                Failure = [Failure jjjj];
                Params(jjjj,:) = -1/10;
            else
                ind_thres = int32(temp + Ts*window -1);
                peak = min(temp1(jjjj,ind_thres:int32(Ts*2*window)));
                ind_peak = find(temp1(jjjj,ind_thres:int32(Ts*2*window)) == peak,1,'first') + ind_thres-1;
            end
        end
        if (~isempty(temp))
            status = 0;
           peak  = abs(peak - mean(mean(temp1(1:int32(ind_thres -Ts*1e-3))))); %% subtract baseline from peak

            %% pchip from 0.5ms before threshold to peak
            PSPface = temp1(jjjj,int32(ind_thres -Ts*0.5e-3):ind_peak);
            step = .1; xs = [1:step:size(PSPface,2)]';
            p = pchip([1:size(PSPface,2)],PSPface,xs);
            halfpeak = find(p>(peak/2),1,'first')*step;
            if ~isempty(halfpeak)
                ind_HPK = halfpeak + int32(ind_thres -Ts*0.5e-3) -1; %% find half peak (slope should be high there so it may be a good place for measuring jitter
            else
                status = -2/10;
                ind_HPK = 0;
            end
            %% Catch stuff that is clearly not clean monosynaptic
            %% Find LOCAL peak
            temp2 = diff(temp1(jjjj,int32(ind_thres -Ts*0.35e-3):ind_peak-int32(0.25e-3*Ts)));
             if temp2(1) <0
                if find(temp2>=0)
                    status = -2/10; %% LOCAL PEAK FOUND
                end
            else
                if find(temp2 <= 0)
                    status =-2/10; %% LOCAL PEAK FOUND
                end
            end
            
            switch (status)
                case 0
                    Params(jjjj,1) = peak;
                    Params(jjjj,2) = ind_HPK;
                case -2/10
                    Params(jjjj,1) = status;
                    Params(jjjj,2) = status;
            end

        end
    catch

    end


    if bdebug
        figure(99);
        hold off
%         plot(intraData1(1:Ts*25e-3,jjjj),'b'); %% this line doesn't work if ~bIntra
        plot(1/Ts*1000.*[1:size(intraData1(1:Ts*2*window,jjjj),1)]-10,intraData1(1:Ts*2*window,jjjj),'b'); %% this line doesn't work if ~bIntra
        hold on
        plot(1/Ts*1000.*[1:size(temp1(jjjj,1:Ts*2*window),2)]-10,temp1(jjjj,1:Ts*2*window),'r');
        line([ind_thres ind_thres]*1/Ts*1000,[min(temp1(jjjj,:))  max(temp1(jjjj,:))],'Color','k')
        line([ind_peak ind_peak]*1/Ts*1000,[min(temp1(jjjj,:))  max(temp1(jjjj,:))],'Color','g')
        line([ind_HPK ind_HPK]*1/Ts*1000,[min(temp1(jjjj,:))  max(temp1(jjjj,:))],'Color','r')
        title([num2str( Params(jjjj,1) ) ' ' num2str(Params(jjjj,2)) ] )
        %         pause;
    end
end;

% %% it is hard to do PSP detection in the presence of spontaneous activity
% 1) baseline fluctuates... if baseline is changing then std will be high and so threshold crossing will
%     probably not occur.  1st derivative threshold could be used to detect change PSPs in this condition
% 2) the max function is used to find the PSP peak. This finds a global maximum, there may be a local maximum that is actually
% from the cell.  It is hard to know which maximum to pick though.  Perhaps one could use the mean latency as a guideline, 
% but you would then have to make assumptions about jitter.


% print('-dpdf', 'tempbass')
 
Params = (Params-(double(Ts*window) -1)) .*dt.*1000;  %% in ms
