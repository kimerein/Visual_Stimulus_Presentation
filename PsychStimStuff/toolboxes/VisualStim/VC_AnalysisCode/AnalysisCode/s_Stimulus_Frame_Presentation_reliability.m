% test for dropped frames.
% STF(1).filename = 'C:\Documents and Settings\Bass\My Documents\Scanziani Lab\VC proj\Data\test frame timing\Data2009-07-09M1_35'; 
STF(1).filename = 'C:\Documents and Settings\Bass\My Documents\Scanziani Lab\VC proj\Data2009-07-17M1_011'; 

for m = 1:length(STF)
    LOADFILE = [STF(m).filename];
       [dataread time]  = daqread([LOADFILE '.daq'], 'Channels', 20, 'Dataformat','native');
       dt = time(2)-time(1);
       
end

temp = diff(dataread>max(dataread)/2);
ind_fall = time(find(temp==-1));
ind_rise = time(find(temp==1));

NFRAMES = length(ind_fall);


temp = diff(ind_rise);
figure(1);  set(gcf,'Position',[356         468        1349         451])
subplot(1,4,1);plot(temp);
title([num2str(NFRAMES) '  ' num2str(time(end)) 's'])
xlabel('frames')

subplot(1,4,2);hist(temp*1000,20)
s= sprintf('\t\t mean min max\ncode+frame jitter: %s ms', num2str([mean(temp) min(temp) max(temp)]*1000));
title(s);
xlabel('ms')

subplot(1,4,4);
plot(time,dataread);
xlabel('s')

temp = ind_fall - ind_rise;
subplot(1,4,3);
hist(temp*1000,20);
s = sprintf('\t\t mean min max\nframe time jitter: %s ms', num2str([mean(temp) min(temp) max(temp)]*1000))
title(s);
xlabel('ms')




%%

SPEC_NFRAMES =  1800 % Specified number of frames
SPEC_FRAMERATE = 30 %Hz


