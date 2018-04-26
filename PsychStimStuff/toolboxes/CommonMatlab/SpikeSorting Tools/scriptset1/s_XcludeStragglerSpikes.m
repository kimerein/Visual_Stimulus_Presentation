Xstd = 2;
for i=1:size(unit,2) % Include only spikes within Xstd of mean
    % find distance of each spike from mean
  unit(i).dist=sum(mean(unit(i).spikes,2)*ones(1,size(unit(i).spikes,2))-double(unit(i).spikes),1);
  s = std(unit(i).dist);
  m = mean(unit(i).dist);
  spike_ind = find(abs(unit(i).dist)<abs(m-Xstd*s)); %find spikes within Xstd of mean
end
