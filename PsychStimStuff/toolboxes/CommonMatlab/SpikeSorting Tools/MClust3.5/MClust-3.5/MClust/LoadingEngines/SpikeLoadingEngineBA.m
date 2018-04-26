function [t, wv]=SpikeLoadingEngineBA(filename,records_to_get,flags);
%function [t, wv]=SpikeLoadingEngineBA(filename,records_to_get,flags);
% should meet Mcluyst 3.5 LoadingEngine specs
load(filename);

if nargin==1; 
t=EG.eventTimes;
wv=EG.wave;
end

if exist('flags')
    switch(flags)

        case 1

        case 2
            t=EG.eventTimes(records_to_get)';
            wv=EG.wave(records_to_get,:,:);
        case 3
            t=EG.eventTimes(records_to_get(1):records_to_get(2))';
            wv=EG.wave([records_to_get(1):records_to_get(2)],:,:);
        case 4
            t=EG.eventTimes(records_to_get(1):records_to_get(2))';
            wv=EG.wave([records_to_get(1):records_to_get(2)],:,:);
        case 5
            t=size(EG.wave,1);

        otherwise
            t=EG.eventTimes';
            wv=EG.wave;
    end
end
