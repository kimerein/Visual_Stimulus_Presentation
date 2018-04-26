function [spiketime wv]= xSpikes4Mclust(abfstruct,EEgroup, EEgroupsize)
% % function xSpikes4Mclust(abfstruct,EEgroup, EEgroupsize)
%     xSpikes4Mclust takes data (loaded from abffile) and extracts spikes for each EEgroup. by applying a threshold
%     on each of the electrodes in EEgroup.  (spikes that occur within minTime are rejected)
% Output: 
%     wfmclust: Spike waveform is formated into  [nspikes x EEgroupsize x 32] for Mclust input
%     spiketime: nspikes x 1
% 
% INPUT:
% abfstruct: %as extracted from abf file with importStruct_abf
% EEgroup: define tetrodes/sterotrode as a [1x4] vector where each value is the position in data.(not ADC number)
%         e.g. [2 3 4 0]
% EEgroupsize: chose 1-N (Mclust only supports upto 4 electrodes)

if EEgroupsize~=size(EEgroup,2)
    error('not all electrodes in EEgroup defined: EEgroupsize ~=size(EEgroup,2)')
end
indexoffset=1;
nIntracellsignals =1;
minTime = 2e-3; % 2ms
overlapWin = minTime/abfstruct.dt;
std_threshold =7; %number of std for threshold
maxpeakamp = 10000*std_threshold;  %in std
maxpeakwidth =  2e-3;
index_beg_spike = [];

% for i=1:size(EEgroup,2) %% each EEgroup
temp = [];
for j=1:EEgroupsize  % all Extracellular signals , threshold
    triggerEE = EEgroup(j);
    if triggerEE >0

        if ~isinteger(abfstruct.data)
            fstd(triggerEE) =  mean(std(abfstruct.data(:,indexoffset+triggerEE:abfstruct.Nchan:end)));
        else
            fstd(triggerEE) =  mean(int16std(abfstruct.data(:,indexoffset+triggerEE:abfstruct.Nchan:end)));
        end

        fmean(triggerEE) =  mean(mean(abfstruct.data(:,indexoffset+triggerEE:abfstruct.Nchan:end)));

        %%% CHANGE code so spikes are in std units
        %% REmember threshold is mean - std*std_threshold
        tic
        temp = [temp; int32(FindSpikes2(abfstruct.data(:,indexoffset+triggerEE:abfstruct.Nchan:end), (fmean(triggerEE)-(std_threshold*fstd(triggerEE))) ,1/abfstruct.dt,0))];
        toc
        %     nSpikes(triggerEE) = size(temp,1);
        index_beg_spike = [index_beg_spike; temp];
        %     stemp = sprintf('N: %d Spike Detection Parameters\tSTD:\t%1.3f, M:\t%1.3f,\tTH:%d',nSpikes(triggerEE), fstd(triggerEE),fmean(triggerEE),std_threshold)
    end
end


if EEgroupsize > 1 %% define neighbor electrodes for xSpikesEE2
    neighbors = [EEgroup(2:end)]
else
    neighbors = 0;
end
neighbors(find(neighbors ==0)) = neighbors(find(neighbors ~=0,1,'first'));


% spikes must occur > then overlapWin (# of samples) to be accepted
index_beg_spike = sort(index_beg_spike);
index_beg_spike = index_beg_spike([1 (find(diff(index_beg_spike)>overlapWin)'+1)]);
[allElectrodesXspikes spikeData temp skipped]  = xSpikesEEmclus(index_beg_spike,abfstruct.data,EEgroup(1),1/abfstruct.dt,abfstruct.Nchan,maxpeakamp+fmean(EEgroup(1)),maxpeakwidth,neighbors,nIntracellsignals,fstd);
%% zero electrodes that are zero
% spiketime = spikeData(find(spikeData(:,1)~=0),1)';
spiketime = spikeData(:,1);

 wv = reshape(allElectrodesXspikes,32,EEgroupsize,size(allElectrodesXspikes,2));
 wv = permute(wv,[3 2 1]);
 if find(EEgroup==0)
     wv(:,find(EEgroup==0),:) = 0;
 end

