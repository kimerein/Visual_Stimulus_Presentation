  
function timeOfFrame = checkFrames(PATH,STF,EFRAMES,chn)
% function timeOfFrame =  checkFrames(PATH,STF,EFRAMES,chn)
% 
% function reads chn with frame info from daq (STF) and plots number
% of high-low transisions (should be frames) an their jitter
% BA080409

if nargin<4
chn = 20; % channel with TTL frame information default
end
timeOfFrame = zeros(length(STF)*EFRAMES,1,'single'); % declare variable

for m = 1:length(STF)
        LOADFILE = [PATH STF(m).filename];
        if STF(1).MAXTRIGGER>1;  error('MaxTrigger must be 1'); % not writtne for more than one trigger per trial TrigRange = [1 STF(m).MAXTRIGGER];
    else TrigRange = STF(1).MAXTRIGGER; end

    [data dt time] = loadDAQData(LOADFILE,chn,TrigRange);
    time = single(time);
    if m == 1 swlength = length(data)*dt;   end% sec
    
    temp = diff(data>max(data)/2);
    ind_fall = time(find(temp==-1));
    NFRAMES = length(ind_fall); 

    if 1 % plot jitter of frames
        ind_rise = time(find(temp==1));
        
        temp = diff(ind_fall);
        figure;  set(gcf,'Position',[ 1681         582         566         402])
        subplot(1,3,1);plot(temp);
        title(['N_{fm}:' num2str(NFRAMES) '  T_{total} ' num2str(time(end)) 's'])
        xlabel('frames'); axis tight;
        ylabel('fm-fm time(s)'); axis tight;
        
        subplot(1,3,2);hist(temp*1000,20)
        s= sprintf('HGRAM fm jitter\nmean min max\ncode+frame jitter: %s ms', num2str([mean(temp) min(temp) max(temp)]*1000));
        title(s);
        xlabel('ms'); axis tight;
        
        subplot(1,3,3);
        plot(time,data);      xlabel('s'); axis tight;
        title('raw signal');axis tight;
        
        temp2 = ind_fall(find(temp == max(temp))); % plot a red line at largest frame
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

    % concatenate date from conseq files
    timeOfFrame(1+EFRAMES*(m-1):m*(EFRAMES))= [ind_fall+(m-1)*swlength];

    clear data time ind_fall temp;
end