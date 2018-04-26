function [params] = extractWOI(data,ind,WOI)
% function extractWOI(data,ind,frac)
%
% returns WOI around each ind within data. data is treated as a linear
% array
%
% INPUT: data  = data should contain 1 channel data from an experiment
%                each row should contain 1 sweep
%      NOTE: data should be filtered because diff will be used to find
%      peaks and troughs in data
%
%       ind = index 
%       WOI = 1x2 array specifying the window around index
%                     [before after]
% OUTPUT:
%         params(1).ind_fracpeak = index where frac occurs for each ind
%             if WOId does not exist for a ind (index in experiment)
%         params  = struct containing
%               params(1).ind_insw =  index in sweep
%               params(1).risetime = peak to trough
%               params(1).eventsw = sweep number; %% MEM make into int16
%               params(1).max2min = amplitude (pk to pk)
%               params(1).skipped  = 1 if is skipped
% BA103006
bdebug =0;

ddata = diff(data);

Ncycles = length(ind);


%% DECLARE PARAMS
params = struct([])
params(1).ind_insw = zeros(1,Ncycles); %% MEM make into int32
params(1).risetime = params(1).ind_insw ;
params(1).eventsw = zeros(1,Ncycles); %% MEM make into int16
params(1).max2min = zeros(1,Ncycles);
params(1).ind_fracpeak = params(1).ind_insw;

params(1).skipped = zeros(1,Ncycles,'int16')
ii = 1;
for i = 1: Ncycles
    try
        %         check same sweep
        if floor((ind(i) - WOId(1))/size(data,1)) == floor((ind(i) + WOId(2))/size(data,1))
            %% find where the trough begins by finding were the derivative goes positive
            trough_beg = find(ddata(ind(i)-WOId(1):ind(i)) >0,1,'last') + ind(i)-WOId(1); %% add 1 because of diff
            %% find trough
            if ~isempty(trough_beg) %% no rise trough is found in dWOI.. rise time is too slow
                trough_max = find(ddata(ind(i):ind(i)+round(WOId(2))) > 0,1,'first') + ind(i); %% add 1 because of diff
                if ~isempty(trough_max) %% no trough is found in dWOI.. rise time is too slow
                    %% find amp of frac of peak
                    side = data(trough_beg)-frac*abs(data(trough_max)-data(trough_beg)); %% ATTN: signal pust be all POSITIVE
                    temp1 = min(abs(data(trough_beg:trough_max)-side));

                    params(1).ind_fracpeak(ii) = find(abs(data(trough_beg:trough_max)-side)==temp1(1),1) + trough_beg-1;
                    params(1).eventsw(ii) = floor((params(1).ind_fracpeak(ii) - WOId(1))/size(data,1));
                    params(1).ind_insw(ii) = params(1).ind_fracpeak(ii) - double(params(1).eventsw(ii))*size(data,1);
                    params(1).risetime(ii) = trough_max-trough_beg; %% time from peak before trough to trough
                    params(1).max2min(ii) = abs(data(trough_max)-data(trough_beg)); %% amplitude of cycle from maximum before osc to trough
                    %% waveforms
                    %                  events(ii,:) = data(params(1).ind_fracpeak(ii)-WOI(1):params(1).ind_fracpeak(ii)+WOI(2));
                    %                  eventtrigger(ii,:) = intradata(params(1).ind_fracpeak(ii)-WOI(1):params(1).ind_fracpeak(ii)+WOI(2));



                    if bdebug
                        figure(7)
                        clf
                        plot(data(ind(i)-WOId(1):ind(i)+WOId(2)))
                        %                         plot(data(max(1,ind(i)-WOId(1)*2):ind(i)+WOId(2)*2))
                        %                 plot(abs(data(ind(i)-trough_beg:ind(i)+trough_max)-side))
                        hold on
                        %                 plot((ddata(ind(i)-WOId:ind(i))-mean(mean(ddata,1)))*100)

                        line([trough_beg-(ind(i)-WOId(1)) trough_beg-(ind(i)-WOId(1))], [min(data(ind(i)-WOId(1):ind(i)+WOId(2))) max(data(ind(i)-WOId(1):ind(i)+WOId(2)))])
                        line([trough_max-(ind(i)-WOId(1)) trough_max-(ind(i)-WOId(1))], [min(data(ind(i)-WOId(1):ind(i)+WOId(2))) max(data(ind(i)-WOId(1):ind(i)+WOId(2)))],'Color','r')
                        line([params(1).ind_fracpeak(ii)-(ind(i)-WOId(1)) params(1).ind_fracpeak(ii)-(ind(i)-WOId(1))], [min(data(ind(i)-WOId(1):ind(i)+WOId(2))) max(data(ind(i)-WOId(1):ind(i)+WOId(2)))],'Color','g')
                        title(num2str(params(1).risetime(ii),'%d'));
                        ii
%                         pause;
                    end

                    ii=ii+1;
                else
                    params(1).skipped(i) = 1;
                end
            else
                params(1).skipped(i) = 1;
            end
        else
            params(1).skipped(i) = 1;
        end

        if mod(i,50) ==0 %% Print progress
            sprintf([num2str(i) ' of ' num2str(Ncycles)]);
        end
    catch
        sprintf(['ERROR at i=' num2str(i,'%d') ' of ' num2str(Ncycles)])
        keyboard
    end
end

%% get rid of zeros from skipped
s = sum(params(1).skipped);
params(1).ind_insw = params(1).ind_insw(1:end-s); %% MEM make into int32
params(1).risetime = params(1).risetime(1:end-s) ;
params(1).eventsw = params(1).eventsw(1:end-s) ; %% MEM make into int16
params(1).max2min = params(1).max2min(1:end-s);
params(1).ind_fracpeak = params(1).ind_fracpeak(1:end-s);
