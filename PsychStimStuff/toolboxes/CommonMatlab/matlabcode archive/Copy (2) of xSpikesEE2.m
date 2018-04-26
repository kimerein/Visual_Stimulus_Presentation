function [Xspikes spikeData sumspiketimeSw skipped] = xSpikesEE2(index_beg_spike,data,Chn_ind,Ts,N_chn,maxpeakamp,maxpeakwidth,neighbors_ind,nIntracellsignals,WOI)
% Differs from xSpikesEE in that data is all electrode data here
% in xSpikesEE it is only the data from one ADC.
%INPUT:  data
% and the index of the beginning of putative spikes (extracted with
% FindSpike).
%       Chn_ind = ADC index (starts @ 1)
%
% OUTPUT:
%  Xspikes - extracted spikes Extra
%  spikeData - (index_beg_spike, sweepnum, timeInSweep, InterspikeInterval)
%  sumspiketimeSw -
%       Each column is one sample in time. Value is the number of spikes
%       occuring at that time
%
tic
indexoffset = 1  ; % data starts with time column and Intracellular column

if nargin < 8
    neighbors_ind = [Chn_ind Chn_ind]; % left for backward compatibility
end

if nargin < 9
    nIntracellsignals = 1;
    warning('nIntracellsignals unspecified in xSpikesEE2.m, [nIntracellsignals = 1]');
end
if neighbors_ind == -1 % all electrodes are neigbors
   neighbors_ind = [(nIntracellsignals+1):(Chn_ind-1) (Chn_ind+1):(N_chn)];
end
%PEAK extaction PARAMETERS
%Window of interest around Beg_spike -500us + 2500us
if nargin < 10
    WOI = [1000e-6 1000e-6; 2000e-6 2000e-6]; % left for backward compatibility
end
Waveform_size = [floor((WOI(1)+WOI(2))*Ts+1); floor((WOI(3)+WOI(4))*Ts+1)];

if isinteger(index_beg_spike)
    index_beg_spike = single(index_beg_spike);
end

N_samples = size(data(:,indexoffset+Chn_ind:N_chn:end),1)*N_chn;  % N samples for all channels in 1 sweep
% correct index_beg_spike so it is correct for whole dataset
indexOffsetElectrode = (indexoffset+(Chn_ind-1))*N_samples/N_chn;
for i = 1: size(neighbors_ind,2)
    indexNeighbor(i) = (indexoffset+(neighbors_ind(i)-1))*N_samples/N_chn;
end

% extract electrode data from all data (this step is memory intensive and
% can be avoided if data is used instead.
lastsweepnum =0;
lastPeak1 = -10000;
if ~isempty(index_beg_spike)
    sumspiketimeSw = zeros(1,N_samples/N_chn);
    Nextactedspikes = 0;

    for (ii  = 1:size(index_beg_spike,1)) % for all spikes found
        z = 1;
        sweepnum = floor(index_beg_spike(ii,1)/(N_samples/N_chn))+1;
        temp = floor((index_beg_spike(ii,1)+ (WOI(2)+1e-3)*Ts)/(N_samples/N_chn))+1;
        temp2 = floor((index_beg_spike(ii,1)- (WOI(1))*Ts)/(N_samples/N_chn))+1;

        indexb = (sweepnum-1) * (N_samples/N_chn)*(N_chn-1);
        electrode_index_beg_spike = indexOffsetElectrode + indexb + index_beg_spike(ii,1);

        if ((temp == sweepnum & temp2 == sweepnum)) %check window around index_beg_spike exists in the same sweep
try
            temp2 = electrode_index_beg_spike+ 1e-3*Ts;  % min in spike must occur within 1ms of crossing threshold
            temp1 = data (electrode_index_beg_spike:temp2)';
catch %DEBUG
    a = ii
end
            if( min(temp1)> maxpeakamp) %
                if false %DEBUG
                    plot(temp1);
                    pause;
                end
                temp1 = find(temp1 == min(temp1));
                Peak1 = temp1(1,1) + electrode_index_beg_spike;
                if( Peak1 - lastPeak1 )> Ts*1e-3   % must be at least 1ms since last spike

                    for i=1: size(neighbors_ind,2)
                        neighbor_index_peak_spike_index(i) = indexNeighbor(i) + indexb + index_beg_spike(ii,1) + temp1(1,1);
                    end

                    spiketime = temp1(1,1) + index_beg_spike(ii,1) ;
                    timeInsweep = spiketime -(sweepnum-1)*N_samples/N_chn;

                    if (timeInsweep > 0 ) %% (what is this condition?
                        Nextactedspikes = Nextactedspikes +1;

                        %WAVEFORM
                        % assume 2 neighbors (otherwise must do for loop that
                        % copies the data and would be slower)
                        try  %for DEBUGGING easily
                                                    if ii == 5438% DEBUG
                                                        a = ii
                                                    end
                                                    % MUST BE MANUALLY CHANGED for more Neighbors
                            Xspikes(:,Nextactedspikes,z) = [full(data(Peak1(1,1)-WOI(1)* Ts:Peak1(1,1)+WOI(2)* Ts)) ...
                                data(neighbor_index_peak_spike_index(1)-WOI(1)* Ts:neighbor_index_peak_spike_index(1)+WOI(2)* Ts) ...
                                data(neighbor_index_peak_spike_index(2)-WOI(1)* Ts:neighbor_index_peak_spike_index(2)+WOI(2)* Ts)];
                            skipped(ii) =0;
                            lastPeak1 = Peak1;
                            if (mod(ii,1000) == 0)
                                Nextactedspikes
                                figure(1);
                                plot(Xspikes(1:int16(Waveform_size(1)*(size(neighbors_ind,2)+1)),Nextactedspikes,z));
                                hold on
                            end

                        catch
                            temp = sprintf('Spike number: %d',ii);
                            error(temp);
                        end

                        if false %DEBUG
                            figure(2);
                            plot(Xspikes(1:int16(round(Waveform_size(1)*(size(neighbors_ind,2)+1))),Nextactedspikes,z));
                            pause;
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

            else
                skipped(ii) = -3;
            end%% Window around spike %%Max amp
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
    stemp = sprintf('EE:\t%d\nNeighbors:\t%d %d\nSpikesX:\t%d of %d \tSkipped: %d,%d,%d,%d\n\t%1.3f sec',Chn_ind,neighbors_ind(1),neighbors_ind(2),...
        Nextactedspikes,size(index_beg_spike,1),Spikeskipped,temp1)