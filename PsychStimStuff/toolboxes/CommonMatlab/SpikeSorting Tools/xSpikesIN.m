function [Xspikes spikeData skipped misalign_ind] = xSpikesIN(index_beg_spike,data,Ts,N_chn)
% [Xspikes spikeData skipped misalign_ind] = xSpikesIN(index_beg_spike,data,Ts,N_chn)
% INPUT:  data (and optionally extracellular data if N_chan ==2)
% and the index of the beginning of putative spikes (extracted with
% FindSpike).
%
% If N_chan >=2 channel 1 is assumed to be intra and channel 2 extra
% all other channels are ignored
%
% OUTPUT:
%  Xspikes - extracted spikes Intra (Extra if extra provided)
%  spikeData - (index_beg_spike, sweepnum, timeInSweep, InterspikeInterval)
%  sumspiketimeSw - 
%       Each column is one sample in time. Value is the number of spikes
%       occuring at that time
%  skipped - holds a code describing why a particular index_beg_spike was
%       no extracted as a spike
%  misalign_ind - holds difference between Peak1(Intra) and Trough1(Extra)

        %PEAK extaction PARAMETERS
        %Window of interest around Beg_spike -500us + 2500us
WOI = [1000e-6 1000e-6; 2000e-6 2000e-6];
Waveform_size = [(WOI(1)+WOI(2))*Ts+1; (WOI(3)+WOI(4))*Ts+1];
maxpeakwidth = 2e-3;

lastsweepnum =0;
N_samples = size(data,1);
if (N_chn >=2)
    intraCell_data = data(:,1+1:N_chn:end);
    extraCell_data = data(:,1+2:N_chn:end);
else
    intraCell_data = data;
end
Xspikes = [];
spikeData = [];
skipped = [];
misalign_ind = [];

if ~isempty(index_beg_spike)
    %Note aligned temporally by 1st zero crossing of Intracellular spike
    sumspiketimeSw = zeros(1,N_samples);
    Nextactedspikes = 0;

    for (ii  = 1:size(index_beg_spike,1)) % for all spikes found
        z = 1;
        sweepnum = floor(index_beg_spike(ii,1)/(N_samples)) +1;
        temp = floor((index_beg_spike(ii,1)+ maxpeakwidth*Ts)/(N_samples)) +1;

        if (temp == sweepnum) %check window around index_beg_spike exists in the same sweep
            try 
                temp1 = intraCell_data (index_beg_spike(ii,1):index_beg_spike(ii,1)+ 1.6e-3*Ts)';
%             if( min(temp1)> maxpeakamp) %
                Peak1 = min(find(temp1 == max(temp1)) + index_beg_spike(ii,1));
            catch
                ii 
            end
                %take first Peak
                % EDGE of SWEEP CHECK
                %check whole waveform is in sweep
                if (mod(Peak1(1,1),size(intraCell_data,1))/Ts > WOI(1)...
                        & (size(intraCell_data,1)- mod(Peak1(1,1),size(intraCell_data,1)))/Ts> WOI(2)) %check if WOI exists around spike (if to early or late in sweep it may not
                    Nextactedspikes = Nextactedspikes +1;

                         %SPIKETIME DATA for ANALYSIS
                    spiketime(Nextactedspikes,1) = index_beg_spike(ii,1);
                    timeInsweep = floor(spiketime(Nextactedspikes,1) -(sweepnum-1)*N_samples);
                    %FIND ISI
                    if (lastsweepnum == sweepnum)
                        isi = timeInsweep - lastsweeptime;
                        lastsweeptime = timeInsweep;
                    else  %first spike of a new sweep
                        lastsweepnum = sweepnum;
                        lastsweeptime = timeInsweep;
                        isi= NaN;
                    end

                    spikeData(Nextactedspikes,:) = [index_beg_spike(ii,1) sweepnum timeInsweep isi];

                    %All sweeps
%                     sumspiketimeSw(1,timeInsweep) = sumspiketimeSw(1,timeInsweep) +1; %Each column is one sample in time. Value is the number of spikes occuring at that time
                    sumspiketimeSw = [];
                    %EXTRACELLULAR SPIKE
                    z = 3;
                    temp1 = extraCell_data (index_beg_spike(ii,1):index_beg_spike(ii,1)+ 1.6e-3*Ts)';
                    Trough1 = min(find(temp1 == min(temp1)) + index_beg_spike(ii,1));
                    Xspikes(1:int16(Waveform_size(1)),Nextactedspikes,z) = extraCell_data (Peak1-WOI(1)* Ts:Peak1+WOI(2)* Ts)';
                    %INTRACELLULAR WAVEFORM
                    
                    z = 1;
                    Xspikes(1:int16(Waveform_size(1)),Nextactedspikes,z) = intraCell_data(Peak1-WOI(1)* Ts:Peak1+WOI(2)* Ts)';
                    
                    
                    misalign_ind(ii) = Peak1 - Trough1;
         
                   
                    if false %DEBUG
                        hold off
                        figure(1)
                        subplot(1,2,1)
                        plot(temp1);
                        hold on
                        subplot(1,2,2)
                        plot(Xspikes(1:int16(Waveform_size(1)),Nextactedspikes,z));

                        pause;
                    end

                else
                    skipped(ii) = -2;
                end%% Window around spike
%             end %%Max amp
        else
            skipped(ii)= -1;
        end %% Centering window for search for PEAK exists
    end % for each spike
end%ISEMPTY BEG_SPIKE
