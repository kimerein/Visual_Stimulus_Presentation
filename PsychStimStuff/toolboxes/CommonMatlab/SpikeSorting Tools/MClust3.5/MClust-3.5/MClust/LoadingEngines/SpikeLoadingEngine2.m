function [t, wv]=SpikeLoadingEngine2(filename,records_to_get,flags);
%function [t, wv]=SpikeLoadingEngine2(filename,records_to_get,flags);
% should meet Mcluyst 3.5 LoadingEngine specs
wv = [];
load(filename);

if nargin==1; 
t=ind_xth;
wv=spikes;
end

if exist('flags')
switch(flags)
    
    case 1
        
    case 2
        t=ind_xth(records_to_get);
        wv=spikes(records_to_get,:,:);
    case 3
        a=find(ind_xth>=records_to_get(1), 1, 'first');
        b=find(ind_xth<=records_to_get(2), 1, 'first');
        t=ind_xth(a:b);
        wv=spikes([a:b],:,:);
    case 4
        t=ind_xth(records_to_get(1):records_to_get(2));
        wv=spikes([records_to_get(1):records_to_get(2)],:,:);
    case 5
        t=size(spikes,1);
        
    otherwise
        t=ind_xth;
wv=spikes;
end
end

if ~isempty(wv)
    for i=1:size(wv,2)
        m(i) = abs(min(min(wv(:,i,:))));
    end
    wv =wv ./ repmat(m,[size(wv,1),1,size(wv,3)]);
end