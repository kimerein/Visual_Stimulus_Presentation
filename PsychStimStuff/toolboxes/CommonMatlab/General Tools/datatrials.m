function dataout = datatrials(data,swl)
% function dataout = datatrials(data,swl)

sw = floor(data./swl) ; % find which sweep a spike occurs in
data = data - sw*swl; % zero so that each spiketime is in reference to the correct episode (not beginning of the experiment)
Nsw = unique(sw);
for i=1:length(Nsw)
    dataout(i).times = data(find(sw==Nsw(i)));
end