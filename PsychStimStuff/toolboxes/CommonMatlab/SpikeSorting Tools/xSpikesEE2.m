function [Xspikes spikeData sumspiketimeSw skipped] = xSpikesEE2(index_beg_spike,data,Chn_ind,Ts,N_chn,maxpeakamp,maxpeakwidth,neighbors_ind,nIntracellsignals,WOI,mean, std)
% Differs from xSpikesEE2 in that data is all electrode data here
% in xSpikesEE it is only the data from one ADC.
% xSpikeEE2 works for intracellular spikes as well.  The Intracellular
% spikes are aligned to peaks rather then troughs this happens if Chn_ind
% (the input channel in ADC 1 based units) is <= nIntracellsignals

%INPUT:  data
% and the index of the beginning of putative spikes (extracted with
% FindSpike).
%       Chn_ind = ADC index (starts @ 1)
%       std: vector of std for each signal (used to determine which signal
%       data should be aligned by.
%
% OUTPUT:
%  Xspikes - extracted spikes Extra
%  spikeData - (index_beg_spike, sweepnum, timeInSweep, InterspikeInterval)
%  sumspiketimeSw -
%       Each column is one sample in time. Value is the number of spikes
%       occuring at that time
%
tic
indexoffset = single(1)  ; % data starts with time column and Intracellular column
bAllEE = 0;
bdebug = 0;

Nextactedspikes = 0;
if nargin < 8
    n_neighbors = 0;
end

if nargin < 9
    nIntracellsignals = 1;
    warning('nIntracellsignals unspecified in xSpikesEE2.m, [nIntracellsignals = 1]');
end

n_neighbors = size(neighbors_ind,2);
if neighbors_ind == -1 % all electrodes are neigbors
    neighbors_ind = [(nIntracellsignals+1):(Chn_ind-1) (Chn_ind+1):(N_chn)]; %% ADC number starting at 1
    bAllEE = 1;
    n_neighbors = single(size(neighbors_ind,2));
elseif neighbors_ind ==0
    n_neighbors = 0;
    neighbors_ind = [];
end
%PEAK extaction PARAMETERS
%Window of interest around Beg_spike -500us + 2500us
if nargin < 10
    WOI = double([750e-6 750e-6; 2000e-6 2000e-6]); % left for backward compatibility
end
WOI= double(WOI); N_chn = double(N_chn); Ts = double(Ts);

Waveform_size = [floor((WOI(1)+WOI(2))*Ts+1); floor((WOI(3)+WOI(4))*Ts+1)];

if isinteger(index_beg_spike)
    index_beg_spike = single(index_beg_spike);
end

N_samples = double(size(data(:,indexoffset+Chn_ind:N_chn:end),1)*N_chn);  % N samples for all channels in 1 sweep
% correct index_beg_spike so it is correct for whole dataset
indexOffsetElectrode = double((indexoffset+(Chn_ind-1))*N_samples/N_chn);
for i = 1: n_neighbors
    indexNeighbor(i) = double((indexoffset+(neighbors_ind(i)-1))*N_samples/N_chn);
end

% extract electrode data from all data (this step is memory intensive and
% can be avoided if data is used instead.
lastsweepnum =0;
lastPeak1 = -10000;

%Preallocate space for variables (speeds up code)
if isinteger(data)
    Xspikes = zeros((int32(WOI(1)*Ts)+int32(WOI(2)*Ts)+1)*(n_neighbors+1),size(index_beg_spike,1),'int16');
else
    Xspikes = zeros((int32(WOI(1)*Ts)+int32(WOI(2)*Ts)+1)*(n_neighbors+1),size(index_beg_spike,1),'single');
end
skipped = int8(zeros(1,size(index_beg_spike,1)));
spikeData = int32(zeros(size(index_beg_spike,1),4));
sumspiketimeSw = zeros(1,N_samples/N_chn,'int32');

if ~isempty(index_beg_spike)
    Nextactedspikes = 0;

    for (ii  = 1:size(index_beg_spike,1)) % for all spikes found
        z = 1;
        sweepnum = double(floor(index_beg_spike(ii,1)/(N_samples/N_chn))+1);
        temp = double(floor((index_beg_spike(ii,1)+ ((WOI(2)+1e-3)*Ts))/(N_samples/N_chn))+1);
        temp2 = double(floor((index_beg_spike(ii,1)- (WOI(1)*Ts))/(N_samples/N_chn))+1);

        if ((temp == sweepnum & temp2 == sweepnum)) %check window around index_beg_spike exists in the same sweep
            % determin index
            indexb = double(((sweepnum-1)) * (N_samples/N_chn)*(N_chn-1));
            electrode_index_beg_spike = double(indexOffsetElectrode + indexb + index_beg_spike(ii,1));

            %define window for peak search
            PreWindow =  (0e-3);
            temp2 = double(electrode_index_beg_spike+  (1e-3*Ts));  % min in spike must occur within 1ms of crossing threshold
            temp1 = double(data ((electrode_index_beg_spike-PreWindow*Ts):temp2))';
%             Amp = abs(min(temp1)-mean(Chn_ind))/std(Chn_ind);
%             Amp = abs(min(temp1)-mean(Chn_ind));
            align_chn = Chn_ind;
            if false %DEBUG
                plot(temp1);
                pause;
            end
            
            %% ADD spline interpolation to find peak (noise at the peak can lead to missalignment of
            %% peaks)
            if Chn_ind <= nIntracellsignals %% signal is intracellular
                min_index = (find(temp1 == max(temp1),1,'first'));
            else
                min_index = (find(temp1 == min(temp1),1,'first'));
            end
  
            if bAllEE %% align to largest signal if triggers came from all electrodes
                for(w=1:n_neighbors) % find the largest signal (relative to std)
                    neighbor_index_beg_spike_index(w) = double(indexNeighbor(w) + indexb + index_beg_spike(ii,1));
                    wtemp(w) = neighbor_index_beg_spike_index(w)+  1e-3*Ts;  % min in spike must occur within 1ms of crossing threshold
                    temp1 = double(data((neighbor_index_beg_spike_index(w)-double(PreWindow*Ts)): wtemp(w)))';
                    temp_index(w) = min(temp1);
%                     wtemp(w) = abs(temp_index(w)-mean(w))/std(neighbors_ind(w));
                    wtemp(w) = abs(temp_index(w)-mean(w));
                    temp_index(w) = find(temp1 == temp_index(w),1,'first');
                end
                temp3= max(wtemp);
                if temp3 > Amp
                    align_chn = find(wtemp == temp3,1,'first');
                    Amp = temp3;
                    min_index = temp_index(align_chn);
                end
            end
            min_index = double(min_index + PreWindow*Ts);
            %             if((Amp*-1)> maxpeakamp) %

            Peak1 = double(min_index + electrode_index_beg_spike);
            if( Peak1 - lastPeak1 )> Ts*1e-3   % must be at least 1ms since last spike
                for i=1: size(neighbors_ind,2)
                    neighbor_index_peak_spike_index(i) = double(indexNeighbor(i) + indexb + index_beg_spike(ii,1) + min_index);
                end

                spiketime = min_index + index_beg_spike(ii,1) ;
                timeInsweep = spiketime -(sweepnum-1)*N_samples/N_chn;

                if (timeInsweep > 0 ) %% (what is this condition?
                    Nextactedspikes = Nextactedspikes +1;

                    %WAVEFORM
                    try  %for DEBUGGING easily
                        clear temp1;
                        if ii == 3% DEBUG
                            a = ii
                        end
                        temp1 = data(int32((Peak1(1,1)-int32(double(WOI(1)* Ts)))):int32(Peak1(1,1)+int32(double(WOI(2)* Ts))));
                        %% UNTESTED
                        for n  = 1: n_neighbors
                            temp = data(int32(double(neighbor_index_peak_spike_index(n))-int32(double(WOI(1)* Ts))):int32(double(neighbor_index_peak_spike_index(n))+int32(double(WOI(2)* Ts))));
                            temp1 = [temp1 temp];
                        end

                        Xspikes(:,Nextactedspikes) = temp1;
                        skipped(ii) =0;
                        lastPeak1 = Peak1;
%                         lastAmp = Amp;
                        if bdebug %DEBUG

                            if (mod(ii,1000) == 0)
                                Nextactedspikes
                                figure(1);
                                plot(Xspikes(1:int16(Waveform_size(1)*(size(neighbors_ind,2)+1)),Nextactedspikes));
                                hold on
                            end
                        end
                    catch
                        temp = sprintf('Spike number: %d',ii);
                        error(temp);
                    end

                    if bdebug %DEBUG
                        figure(2);
                        plot(Xspikes(:,Nextactedspikes));
                        pause;
                        axis tight;
                    end
                    %SPIKETIME DATA for ANALYSIS
                    %FIND ISI
                    if (lastsweepnum == sweepnum)
                        isi = timeInsweep - lastsweeptime;
                        lastsweeptime = timeInsweep;
                    else  %first spike of a new sweep
                        lastsweepnum = sweepnum;
                        lastsweeptime = timeInsweep;
                        isi= NaN;
                    end

                    spikeData(Nextactedspikes,:) = [spiketime sweepnum timeInsweep isi];
                    % make into struct

                    %All sweeps
                    sumspiketimeSw(1,timeInsweep) = sumspiketimeSw(1,timeInsweep) +1; %Each column is one sample in time. Value is the number of spikes occuring at that time
                end% time in sweep > = 0
            else  % spike occured less then 1ms since last spike
                skipped(ii) = -4;
            end

            %             else
            %                 skipped(ii) = -3;
            %             end%% Window around spike %%Max amp
        else
            skipped(ii)= -1;
        end% %% Centering window for search for PEAK exists
    end % for each spike

end%ISEMPTY BEG_SPIKE
temp1 = toc;
%REPORT
for jj = 1:4
    Spikeskipped(jj) = sum(skipped == -1*jj);
end

%% ADD remove zeros from output data

stemp = sprintf('EE:\t%d\nNeighbors:\t%s\nSpikesX:\t%d of %d \tSkipped: %d,%d,%d,%d\n\t%1.3f sec',Chn_ind,num2str(neighbors_ind),Nextactedspikes,size(index_beg_spike,1),Spikeskipped,temp1)