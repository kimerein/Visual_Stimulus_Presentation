function [Xspikes spikeData sumspiketimeSw skipped] = xSpikesEE(index_beg_spike,data,Ts,N_chn,maxpeakamp,maxpeakwidth,WOI)
% INPUT:  data
% and the index of the beginning of putative spikes (extracted with
% FindSpike).
%
% If N_chn > 1 channel 1 is assumed to be extra
%
% OUTPUT:
%  Xspikes - extracted spikes Extra
%  spikeData - (index_beg_spike, sweepnum, timeInSweep, InterspikeInterval)
%  sumspiketimeSw -
%       Each column is one sample in time. Value is the number of spikes
%       occuring at that time
%
%PEAK extaction PARAMETERS
%Window of interest around Beg_spike -500us + 2500us
if nargin < 7
    WOI = [1000e-6 1000e-6; 2000e-6 2000e-6];
end
    Waveform_size = [(WOI(1)+WOI(2))*Ts+1; (WOI(3)+WOI(4))*Ts+1];

N_samples = size(data,1)*N_chn;
lastsweepnum =0;
if ~isempty(index_beg_spike)
     sumspiketimeSw = zeros(1,N_samples/N_chn);
    Nextactedspikes = 0;

    for (ii  = 1:size(index_beg_spike,1)) % for all spikes found
        z = 1;
        sweepnum = floor(index_beg_spike(ii,1)/(N_samples/N_chn)) +1;
        temp = floor((index_beg_spike(ii,1)+ maxpeakwidth*Ts)/(N_samples/N_chn)) +1;

        if (temp == sweepnum) %check window around index_beg_spike exists in the same sweep
             temp2 = index_beg_spike(ii,1)+ int32(1.6e-3*Ts);
            if(temp2 < size(data,1)*size(data,2))
                temp1 = data (index_beg_spike(ii,1):temp2)';
                if( min(temp1)> maxpeakamp) %
                    if false %DEBUG
                        figure;
                        plot(temp1);
                        pause;
                    end
                    Peak1 = find(temp1 == min(temp1)) + index_beg_spike(ii,1);

                    %take first Peak
                    % EDGE of SWEEP CHECK
                    %check whole waveform is in sweep
                    if (mod(Peak1(1,1),size(data,1))/Ts > WOI(1) & (size(data,1)- mod(Peak1(1,1),size(data,1)))/Ts> WOI(2)) %check if WOI exists around spike (if to early or late in sweep it may not
%                         spiketime = index_beg_spike(ii,1);
                        spiketime = Peak1(1,1);
                        timeInsweep = spiketime -(sweepnum-1)*N_samples/N_chn;

                        if (timeInsweep > 0 )
                            Nextactedspikes = Nextactedspikes +1;

                            %WAVEFORM
                            Xspikes(:,Nextactedspikes,z) = full(data(Peak1(1,1)-WOI(1)* Ts:Peak1(1,1)+WOI(2)* Ts)');
                            if false %DEBUG
                                figure;
                                plot(Xspikes(1:int16(Waveform_size(1)),Nextactedspikes,z));
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
                    else
                        skipped(ii) = -2;
                    end%% Window around spike
                end %%Max amp
            end
        else
            skipped(ii)= -1;
        end% %% Centering window for search for PEAK exists
    end % for each spike
end%ISEMPTY BEG_SPIKE
