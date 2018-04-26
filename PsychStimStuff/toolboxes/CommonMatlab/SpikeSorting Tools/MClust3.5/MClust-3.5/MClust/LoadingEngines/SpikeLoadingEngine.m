function [t, wv]=SpikeLoadingEngine(filename,records_to_get,flags);
%function [t, wv]=SpikeLoadingEngine(filename,records_to_get,flags);
% should meet Mcluyst 3.5 LoadingEngine specs
load(filename);

if nargin==1; 
t=LFPev.ind_fracpeak;
wv=spikes;
end


switch(flags)
    
    case 1
        
    case 2
        t=LFPev.ind_fracpeak(records_to_get);
        wv=spikes(records_to_get,:,:);
    case 3
        a=find(LFPev.ind_fracpeak>=records_to_get(1), 1, 'first');
        b=find(LFPev.ind_fracpeak<=records_to_get(2), 1, 'first');
        t=LFPev.ind_fracpeak(a:b);
        wv=spikes([a:b],:,:);
    case 4
        t=LFPev.ind_fracpeak(records_to_get(1):records_to_get(2));
        wv=spikes([records_to_get(1):records_to_get(2)],:,:);
    case 5
        t=size(spikes,1);
        
    otherwise
        t=LFPev.ind_fracpeak;
wv=spikes;
end

t= t'; % make t n x 1 (where n is the number of spikes