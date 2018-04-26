function [Xspikes spikeData sumspiketimeSw skipped] = xSpikesEEgroup(EEgroup,index_beg_spike,abfstruct,nIntracellsignals,WOI)
%INPUT: EEgroup (electrode definition) 
%   index_beg_spike
%   abfstruct
%        OPTIONAL
%   nIntracellsignals
%   WOI
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
bdebug = 1;
Nextactedspikes = 0;

if nargin < 4
    nIntracellsignals = 1;
        warning('nIntracellsignals unspecified in xSpikesEE2.m, [nIntracellsignals = 1]');
end
if nargin < 5
   WOI = double([500e-6 500e-6; 1500e-6 1500e-6]); % left for backward compatibility
end
% if nargin < 6
%     maxpeakamp = 1000;
% end
% if nargin < 7
%     maxpeakwidth = 1e-3
% end
Ts = double(1/abfstruct.dt);
N_chn = double(abfstruct.Nchan);
n_electrodes = size(EEgroup,2);

%find mean and std of EEgroup
for i = 1:n_electrodes
    triggerEE = EEgroup(i);
    if ~isinteger(abfstruct.data)
        fstd(i) =  mean(std(abfstruct.data(:,indexoffset+triggerEE:N_chn:end)));
    else
        fstd(i) =  mean(int16std(abfstruct.data(:,indexoffset+triggerEE:N_chn:end)));
    end
    fmean(i) =  mean(mean(abfstruct.data(:,indexoffset+triggerEE:N_chn:end)));
end


%PEAK extaction PARAMETERS
%Window of interest around Beg_spike -500us + 2500us
Waveform_size = [floor((WOI(1)+WOI(2))*Ts+1); floor((WOI(3)+WOI(4))*Ts+1)];

if isinteger(index_beg_spike)
    index_beg_spike = single(index_beg_spike);
end

N_samples = double(size(abfstruct.data(:,indexoffset+1:N_chn:end),1)*N_chn);  % N samples for all channels in 1 sweep
% correct index_beg_spike so it is correct for whole dataset
for i = 1: n_electrodes
    indexElectrode(i) = double((indexoffset+(EEgroup(i)-1))*N_samples/N_chn);
end

% extract electrode data from all data (this step is memory intensive and
% can be avoided if data is used instead.
lastsweepnum =0;

%Preallocate space for variables (speeds up code)
if isinteger(abfstruct.data)
    Xspikes = zeros((int32(WOI(1)*Ts)+int32(WOI(2)*Ts)+1)*(n_electrodes),size(index_beg_spike,1),'int16');
else
    Xspikes = zeros((int32(WOI(1)*Ts)+int32(WOI(2)*Ts)+1)*(n_electrodes),size(index_beg_spike,1),'single');
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

            indexb = double(((sweepnum-1)) * (N_samples/N_chn)*(N_chn-1));
            PreWindow =  (0e-3);            %define window for peak search
            for(w=1:1) % find the largest signal (relative to std) %% ALIGNED to first electrode because somehow aligning to the largest one doesn't work
                % determin index
                %% ADD spline interpolation to find peak (noise at the peak can lead to missalignment of
                %% peaks)
                electrode_index_beg_spike_index(w) = double(indexElectrode(w) + indexb + index_beg_spike(ii,1));
                wtemp(w) = electrode_index_beg_spike_index(w)+  1e-3*Ts;  % min in spike must occur within 1ms of crossing threshold
                temp1 = double(abfstruct.data((electrode_index_beg_spike_index(w)-double(PreWindow*Ts)): wtemp(w)))';
                temp_index(w) = min(temp1);
                %                     wtemp(w) = abs(temp_index(w)-mean(w))/std(EEgroup(w));
                wtemp(w) = abs((temp_index(w)-fmean(w))/fstd(w));
                temp_index(w) = find(temp1 == temp_index(w),1,'first');
            end
            temp3= max(wtemp);
            align_chn = find(wtemp == temp3,1,'first');
            Amp = temp3;
            min_index =double(temp_index(align_chn)+ PreWindow*Ts);

            if false %DEBUG
                plot(temp1);
                pause;
            end

            %          if((Amp*-1)> maxpeakamp) %
            for i=1: n_electrodes
                electrode_index_peak_spike_index(i) = double(indexElectrode(i) + indexb + index_beg_spike(ii,1) + min_index);
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
                    %% UNTESTED
                    temp1 = [];
                    for n  = 1: n_electrodes
                        temp = abfstruct.data(int32(double(electrode_index_peak_spike_index(n))-int32(double(WOI(1)* Ts))):int32(double(electrode_index_peak_spike_index(n))+int32(double(WOI(2)* Ts))));
                        temp1 = [temp1 temp];
                    end

                    Xspikes(:,Nextactedspikes) = temp1;
                    skipped(ii) =0;
                     lastAmp = Amp;
                    if bdebug %DEBUG

                        if (mod(ii,1000) == 0)
                            Nextactedspikes
                            figure(1);
                            plot(Xspikes(1:int16(Waveform_size(1)*(size(EEgroup,2)+1)),Nextactedspikes));
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
%                     pause;
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

stemp = sprintf('\nEEgroup:\t%s\nSpikesX:\t%d of %d \tSkipped: %d,%d,%d,%d\n\t%1.3f sec',num2str(EEgroup),Nextactedspikes,size(index_beg_spike,1),Spikeskipped,temp1)