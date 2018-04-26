function [stimes sw swstimes conWV meanWV] = concSpikeData(datafile)
% function [stimes sw swstimes conWV meanWV] = concSpikeData(datafile)
% datafile(N) - struct containing spiketimes and Waveform (WV) (optinal)
% for N data file in an experiment Data containing 
%
% 1x2 struct array with fields:
%     filename
%     spiketimes
%     dt
%     maxsample
%     MAXTRIGGER
%     params
%     LOADPATH

swlength = datafile.maxsample;
dt = datafile.dt;
%concatenate time sesimultaneous 
bSUBTRACT_MEAN = 1; % remove common mode noise by subtracting mean of all tetrodewaveforms.
conWV = [];meanWV = [];stimes = [];sw = [];swstimes = []; 
% DO declare all these variables in there full size first.

for m = 1:length(datafile)
    for j = 1:length(datafile(m).spiketimes)
        
        if nargout >=4
            conWV = [conWV reshape(datafile(m).spiketimes(j).WV,size(datafile(m).spiketimes(j).WV,1)*size(datafile(m).spiketimes(j).WV,2),size(datafile(m).spiketimes(j).WV,3))];
            if bSUBTRACT_MEAN % subtract mean waveform
                meanWV = [meanWV repmat(squeeze(mean(datafile(m).spiketimes(j).WV,2)),size(datafile(m).spiketimes(j).WV,2),1)];
            else meanWV = 0; end
        end

        fileOffsetsw = 0;
        for i = m-1:-1:1;            fileOffsetsw = fileOffsetsw+ datafile(i).MAXTRIGGER;        end % make different files continous by adding time offset before spiketime of file
%         stimes = [stimes; datafile(m).spiketimes(j).ind_sw_trough+(datafile(m).spiketimes(j).sw-1)*swlength+1/dt+fileOffsetsw*swlength]; % add 1 S so there is no correlation between events on different sweeps
        stimes = [stimes; datafile(m).spiketimes(j).ind_sw_trough+(datafile(m).spiketimes(j).sw-1)*swlength+fileOffsetsw*swlength]; % 
        sw = [sw; ones(length(datafile(m).spiketimes(j).sw),1).*(datafile(m).spiketimes(j).sw+fileOffsetsw)];
        swstimes = [swstimes; ones(length(datafile(m).spiketimes(j).ind_sw_trough),1).*datafile(m).spiketimes(j).ind_sw_trough];
    end
end