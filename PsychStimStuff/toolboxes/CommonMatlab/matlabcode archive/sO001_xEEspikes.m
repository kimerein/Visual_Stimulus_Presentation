% Extact spikes
% bin sweeps onto some size
% count spikes
close all;
clear all;
total_extractedspikes = 0;
processedFiles = '';
Nspikes_in_File = [];

 readdirheader =  'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\';
% readdirheader =  '\\132.239.158.164\Scanziani Lab\Patch Data\';
% readdirheader =  'e:\My Documents\Scanziani Lab\Patch Data\';
readfileheader = ['s2_2005_10_21_'];
savefileheader = ['s2_2005_10_21_'];
readfilenumber  = 3;
[filename] = createaxonfilename(readfileheader,readfilenumber(1));
pathfilename = strcat(readdirheader,filename);
[sfilename] = creatematlabfilename(savefileheader,readfilenumber(1));
savefilename = strcat(readdirheader,sfilename);

Chan = 1; %ONLY use 1 channel
inTs = 1/0.1e-3 ;%Sample rate 
sweeps = -1;
%EXTRACT DATA
tic
% [my_data N_sweeps N_chan N_samples] = import_abf(pathfilename,sweeps,1/inTs)%;,Chan,subsample);
[my_data N_sweeps N_chn N_samples] = import_abf(pathfilename,sweeps,1/inTs);%;,Chan,subsample);
toc
if sweeps ~=-1
    N_sweepsExtracted = size(sweeps,2);
else
    N_sweepsExtracted = N_sweeps;
end

tic

subsample = 1; %
Ts = inTs/subsample; %subsampled 
%DATA
time_data = my_data(1:subsample:end,1);
data =my_data(1:subsample:end,1+Chan:N_chn:end) ; %Data w/o


clear my_data; %save memory
% pack;

threshold = -10; %10
maxpeakamp = -100;
% minpeakwidth = .2e-3; %threshold to threshold
maxpeakwidth = 2e-3; %sec

% FIND spikes
data_GT = (data < threshold) ;

temp = [diff(data_GT); zeros(1,size(data_GT,2))]; %insert a row because diff is small then the matrix it comes from
% clear data_GT;
index_beg_spike = find(temp == 1);
index_end_spike = find(temp == -1);
clear temp;

%% CHOOSE SPIKES with min and maxpeakwidth
if (index_end_spike(1) < index_beg_spike(1))  % CASE Falling edge occurs before rising edge
    index_end_spike = index_end_spike(2:end);
end
if (size(index_end_spike,1) < size(index_beg_spike,1) ) %CASE Rising edge with no falling edge
    index_beg_spike = index_beg_spike(1:end-1);
end
for (i=1:size(index_end_spike,1))
    if index_beg_spike(i) >  index_end_spike(i)
        error('')
        i
    end
end
peakwidth = index_end_spike - index_beg_spike;
mask = (peakwidth < maxpeakwidth*Ts); %Choose peaks with correct width %Indices that meet criteria have value 1
index_beg_spike = nonzeros(mask .* index_beg_spike);
NextactedSpikes = size(index_beg_spike,1);
clear mask;clear peakwidth;clear index_end_spike;

%% END of CHOOSE SPIKES
toc;

tic
%% EXTRACTDATA (%TODO funciton)
WOI = [1000e-6 1000e-6; 2000e-6 2000e-6];
Waveform_size = [(WOI(1)+WOI(2))*Ts+1; (WOI(3)+WOI(4))*Ts+1];
extracted_spikes = [];
if ~isempty(index_beg_spike)

    % find Spike interval (beg to beg of spike)
    %             index_beg_spike_sh = circshift(index_beg_spike,-1);
    %             index_beg_spike_sh(size(index_beg_spike_sh,1),1) = 0; %remove element shifted to the end of the vector
    %             interspike_interval = index_beg_spike_sh - index_beg_spike; %interval in index
    %             interspike_interval = circshift(interspike_interval,1);
    %             interspike_interval(1,1) = -1;

    %Note aligned temporally by 1st zero crossing of  spike
    Nextactedspikes = 0;
    sweeptimeSw = zeros(size(index_beg_spike,1),2);
    sumspiketimeSw = zeros(1,N_samples/N_chn);
    for (ii  = 1:size(index_beg_spike,1)) % for all spikes found
        z = 1;
        temp1 = data (index_beg_spike(ii,1):index_beg_spike(ii,1)+ 1.6e-3*Ts)';
        if( min(temp1)> maxpeakamp) %
%             figure(999) %DEBUG
%             plot(temp1);
            
            Peak1 = find(temp1 == min(temp1)) + index_beg_spike(ii,1);
            %take first Peak
            % EDGE of SWEEP CHECK
            SampleInSweep = Peak1(1,1) - floor(Peak1(1,1)/N_samples)*N_samples; %check if WOI exists around spike (if to early or late in sweep it may not
            if ( SampleInSweep >= WOI(1)* Ts && SampleInSweep+WOI(2)* Ts < N_samples) %beginning & end of sweep
                Nextactedspikes = Nextactedspikes +1;
                % NOTE: Nextractedspikes <= size(spiketime,2) because some spikes may have
                % been need sweep edges and not extracted

                %WAVEFORM
                tempXspikes(1:int16(Waveform_size(1)),Nextactedspikes,z) = data(Peak1(1,1)-WOI(1)* Ts:Peak1(1,1)+WOI(2)* Ts)';
%            figure(998) %DEBUG
%             plot(tempXspikes);
%             pause;

                %SPIKETIME DATA for ANALYSIS
                spiketime(Nextactedspikes,1) = index_beg_spike(ii,1);
                %SweepNumber and SweepTime      
                sweepnum = floor(spiketime(Nextactedspikes,1)/(N_samples/N_chn)) +1;
                timeInsweep = spiketime(Nextactedspikes,1) -(sweepnum-1)*N_samples/N_chn;
                spiketimeSw(Nextactedspikes,:) = [sweepnum timeInsweep];
                %All sweeps 
                sumspiketimeSw(1,timeInsweep) = sumspiketimeSw(1,timeInsweep) +1;
                
%                 Sweep =ceil(Peak1(1,1)/N_samples);
%                 if (Sweep > N_sweeps) %HACK not sure why it every can happen
%                     Sweep = N_sweeps;
%                 end
%                 Timeofspike = Peak1(1,1) -N_samples*(Sweep -1);
                
            end
        end


    end % for each spike
    %filename must have 19 columns
    filename = [blanks(19-length(filename)) filename];
    processedFiles = [processedFiles; filename];
%     Nspikes_in_File(ii) = size(index_beg_spike,1)
end%ISEMPTY BEG_SPIKE


extracted_spikes = [extracted_spikes tempXspikes];
toc
% SAVE
save(savefilename,...
    'processedFiles','extracted_spikes','Nextactedspikes','spiketime','N_sweeps','N_sweepsExtracted','N_samples','Ts','subsample','-mat');

figID = 100;
%PLOT Sweep  Threshold & Patched Cell
figure(figID)
subplot(1,2,1)
temp = [1:1];% [1:size(data,2)];
plot(time_data,data(:,temp),'b');
hold all;
stemp = sprintf('%s\n SwFreq(Hz)=%3.0f \nThres= %3.1f  Ts(Hz)= %5.0f',processedFiles,1/(N_samples/N_chn/Ts),threshold,Ts );
ylim([min(data(5000:end,1)) max(data(5000:end,1))]);
title(stemp,'Interpreter','none');
xlabel('Time(s)')
ylabel('Voltage')
% threshold
line([0 time_data(end,1)], [threshold threshold],'Color','k'); %(x1 x2), (y1 y2)
subplot(1,2,2)
% %%PLOT SPIKE WAVEFORMS 
    %TODO turn to histogram
    %   plot criteraila
% figure;
for i=1:size(extracted_spikes,2)
    plot([1:size(extracted_spikes,1)]*1000/Ts,extracted_spikes(:,i));
    hold on
end
stemp = sprintf('%s \nSpike # %d\nthres: %3.1f max:%d maxwidth(ms):%1.1f',processedFiles,Nextactedspikes,threshold,maxpeakamp,maxpeakwidth*1000);
title(stemp,'Interpreter','none');
xlim([0 size(extracted_spikes,1)*1000/Ts]);
ylabel('Voltage');
xlabel('Time (ms)');
% 
%%%%%%%%%%%%%%%%%%%%%%55555555


%PLOT Spiketimes (concatanated sweeps)
%Objective 001
figure(figID+1)
plot(spiketime/Ts,[1:size(spiketime,1)]);
stemp = sprintf('%s \nCum Spike # vs Time \nthres: %3.1f',processedFiles,threshold);
title(stemp,'Interpreter','none');
ylabel('Cum Spike')
xlabel('Time (s)')
% Unity line
line([0 max(spiketime)/Ts],[0 size(spiketime,1)],'Color','k','LineWidth',0.5); 


% MARK SWEEPS beginnings
% clear sweeptime;
% sweepnum = [50 104+33 104+93 284];
% sweeptime = sweepnum*90000;
% sweeptime = [sweeptime; ones(1,size(sweeptime,2))]';
% plot(sweeptime(:,1),sweeptime(:,2),'.r');


% Objective 00

% line(Interval/Ts,[max(spiketimeSw(:,1))-1 max(spiketimeSw(:,1))-1],'Color','r','LineWidth',3); %HORIZONTAL
% line([Interval(1)/Ts Interval(1)/Ts],[min(spiketimeSw(:,1))+1 max(spiketimeSw(:,1))-1],'Color','r','LineWidth',0.5); %VERTICAL
% line([Interval(2)/Ts Interval(2)/Ts],[min(spiketimeSw(:,1))+1 max(spiketimeSw(:,1))-1],'Color','r','LineWidth',0.5); %VERTICAL
