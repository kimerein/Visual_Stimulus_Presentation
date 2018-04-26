function indlight = getLightOnTrials(STF,chn,DAQSAVEPATH)
% function indlight = getLightOnTrials(STF,chn,DAQSAVEPATH)
% to seperate sweeps were light was on from off:
% detects sweeps that have an integratal that is greater than mean/2 of the
% mean integral of all sweeps. and calls those on.
% 
figure(99); clf;
for i = 1:length(STF)
[data dt] = loadDAQData(fullfile(DAQSAVEPATH,STF(i).filename),chn);
temp = sum(data(:,:,1));
m = mean(temp)/2;
indlight{i} = temp> m;
subplot(3,1,1)
  [a x] = hist(temp); h = stairs(x,a); hold all;
  line([m m],[0 1]*max(a),'color',get(h,'color'),'linewidth',2)
subplot(3,1,2);
plot(repmat([1:4:size(data,1)]'*dt,1,sum(indlight{i})),data(1:4:end,indlight{i},1)+max(max(max(data(:,indlight{i},1)),max(ylim)))*1.05,'color',get(h,'color')); hold on;
title('ON')
subplot(3,1,3);
plot(repmat([1:4:size(data,1)]'*dt,1,sum(~indlight{i})),data(1:4:end,~indlight{i},1)+max(max(max(data(:,~indlight{i},1)),max(ylim)))*1.05,'color',get(h,'color')); hold on;
title('OFF')
  sl{i} = STF(i).filename;
end
subplot(3,1,1)
legend(sl);