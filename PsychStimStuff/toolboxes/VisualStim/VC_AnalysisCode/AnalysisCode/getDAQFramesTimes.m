function timeOfFrame = getDAQFramesTimes(DAQSAVEPATH,STF,EFRAMES,DAQchn)
% function timeOfFrame =  getFramesTimes(DAQSAVEPATH,STF,EFRAMES,DAQchn)
%
% function reads DAQchn with frame info from daq (STF) and plots number
% of high-low transisions (should be frames) an their jitter
% INPUT:
%       PATH -
%       STF - struct with field "filename:
%       EFRAMES - scalar with expected number of frames
%       DAQchn - DAQ chn with TTL frames
%
% OUTPUT: timeOfFrame - struct array with field frametime containing frame
% transitions
%
% outputs time of frame (defined as when values goes low)
% in struct swdata.frametime
% BA080409

STA_SAVEPATH = [];
s_RIGSPECIFIC_SPIKESORT_CONFIG

SAVEPATH = fullfile(STA_SAVEPATH,'frametimes');
% save extracted frame times
ind = strfind(STF(1).filename,'_'); 
stemp = fullfile(SAVEPATH, [STF(1).filename(1:ind-1) getfileindex(STF)  '_timeOfFrame']);

if ~isempty(dir([stemp '.mat']))
    load(stemp,'timeOfFrame')
else
    timeOfFrame = getFramesTimeshelper(DAQSAVEPATH,STF,EFRAMES,DAQchn);
    save(stemp,'timeOfFrame'); %% TEMP solution
end

end
function swdata = getFramesTimeshelper(DAQSAVEPATH,STF,EFRAMES,DAQchn)
if nargin<4
    DAQchn = 20; % channel with TTL frame information default
end

    printf('Extracting frame times...');

for m = 1:length(STF)
    LOADFILE = fullfile(DAQSAVEPATH, STF(m).filename);
      % not writtne for more than one trigger per file
    [data dt time] = loadDAQData(LOADFILE,DAQchn,1);
    printf('%s...', STF(m).filename);
    time = single(time);
    if m == 1 swlength = length(data)*dt;   end% sec
    
    temp = diff(data>max(data)/2);
    ind_fall = time(find(temp==-1));
    NFRAMES = length(ind_fall);
    
    if 0 % plot jitter of frames
        ind_rise = time(find(temp==1));
        
        temp = diff(ind_fall); % time between flip calls in stimulus presentation loop
        figure;  set(gcf,'Position',[ 1681         582         566         402])
        subplot(1,3,1);plot(temp*1000);
        title(['N_{fm}:' num2str(NFRAMES) '  T_{total} ' num2str(time(end)) 's'])
        xlabel('frames'); axis tight;
        ylabel('fm-fm time(ms)'); axis tight;
        
        subplot(1,3,2);hist(temp*1000,20)
        s= sprintf('HGRAM fm jitter\nmean min max\ncode+frame jitter: %s ms', num2str([mean(temp) min(temp) max(temp)]*1000));
        title(s);
        xlabel('ms'); axis tight;
        
        subplot(1,3,3); % to do place this on same axis as plottemp
        plot(time,data);      xlabel('s'); axis tight;
        title('raw signal');axis tight;
        
        temp2 = ind_fall((temp == max(temp))); % plot a red line at largest frame
        h=line([1 1].*temp2,[min(data) max(data)]);
        set(h,'color','r','linewidth',2);
        
        if length(ind_fall)==length(ind_rise)
            temp = ind_fall - ind_rise;
            subplot(1,4,3);
            hist(temp*1000,20);
            s = sprintf('mean min max\nframe time jitter: %s ms', num2str([mean(temp) min(temp) max(temp)]*1000));
            title(s);
            xlabel('ms')
        end
        
        clear ind_rise;
    end
    
    % Check right number of frames
    if EFRAMES ~= NFRAMES
        warning(['Wrong number of frames. Expected: ' num2str(EFRAMES) 'Actual: ' num2str(NFRAMES)])
    end
    swdata(m).frametime = ind_fall;
    %     % concatenate date from conseq files
    %     timeOfFrame(1+EFRAMES*(m-1):m*(EFRAMES))= [ind_fall+(m-1)*swlength];
    
    clear data time ind_fall temp;
end
end
