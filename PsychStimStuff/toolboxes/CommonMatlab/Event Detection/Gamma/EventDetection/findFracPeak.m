function [params] = findFracPeak(data,ind_xthres,frac,WOId,bdebug)
% function [params] = findFracPeak(data,ind_xthres,frac,WOId,bdebug)
%
% WARNING: because this alg uses the first derivative to find the trough
% and peaks, events that have a double peak (for example polysynaptic
% currents), may be aligned to the trough within between these two peaks
% and not the actual onset of the current. Practically this is a small
% problem in most data. One can see it by looking at the mean event. A bump before
% the onset of the "main" current suggests that some times there is misalignment.
%`
% Finds index where amplitude is frac*Peak - Trough
%       Peak before ind_xthres
%       Trough after ind_xthres.
%
% INPUT: data  = data should contain 1 channel data from an experiment
%                each row should contain 1 sweep
%      NOTE: data should be filtered because diff will be used to find
%      peaks and troughs in data
%
%       ind_xthres = index of threshold crossing (from thresh function
%       frac  = fraction of trough to peak to find (0<frac<=1
%       WOId = 1x2 array specifying the window to look for a trough
%                     after  ind_xthres and for a peak before ind_xthres
%                     [before after]
% OUTPUT:
%         params(1).ind_fracpeak = index where frac occurs for each ind_xthres
%             if WOId does not exist for a ind_xthres (index in experiment)
%         params  = struct containing
%               params(1).ind_insw =  index in sweep
%               params(1).risetime = peak to trough
%               params(1).eventsw = sweep number; %% MEM make into int16
%               params(1).max2min = amplitude (pk to pk)
%               trough
%               params(1).skipped  = 1 if is skipped
% BA103006
%
% BA113001 changed params.risetime &.max2min to be initiated to NaN
%          added Skipped case where trough and peak are detected at the
%             same point.
%          added params.max (amplitude of peak before threshold crossing)
if nargin < 5
    bdebug =0;
end


Ncycles = length(ind_xthres);


%% DECLARE PARAMS
params = struct([]);
params(1).ind_insw = zeros(1,Ncycles); %% MEM make into int32
params(1).risetime = NaN(1,Ncycles) ;
params(1).eventsw = zeros(1,Ncycles); %% MEM make into int16
params(1).max2min = NaN(1,Ncycles);
params(1).ind_fracpeak = params(1).ind_insw;
params(1).max = NaN(1,Ncycles);
params(1).skipped = zeros(1,Ncycles,'int16');
ii = 1;
if isrowvector(data)
    data = data';
end
ddata = diff(data);

for i = 1: Ncycles
    try
        if floor((ind_xthres(i) - WOId(1))/size(data,1)) == floor((ind_xthres(i) + WOId(2))/size(data,1))  %         check same sweep
            %% find where the trough begins by finding were the derivative goes positive
            if ind_xthres(i)> WOId(1) %check first point as WOId(1) before it
                trough_beg = find(ddata(ind_xthres(i)-WOId(1):ind_xthres(i)) >0,1,'last') + ind_xthres(i)-WOId(1); %% add 1 because of diff
                %% find trough
                if bdebug
                    figure(8)
                    plot(data(max(1,ind_xthres(i)-WOId(1)*20):min(ind_xthres(i)+WOId(2)*20,length(data))))

                    figure(7)
                    clf
                    plot([-WOId(1):WOId(2)],data(ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2)))
                    %                     plot(data(ind_xthres(i)-WOId(1)*10:ind_xthres(i)+WOId(2)*10))
                end
                if ~isempty(trough_beg) %% no rise trough is found in dWOI.. rise time is too slow
                    trough_max = find(ddata(ind_xthres(i):ind_xthres(i)+round(WOId(2))) > 0,1,'first') + ind_xthres(i); %% add 1 because of diff
                    if ~isempty(trough_max) %% no trough is found in dWOI.. rise time is too slow
                        %% find amp of frac of peak
                        if trough_beg~= trough_max
                            side = data(trough_beg)-frac*abs(data(trough_max)-data(trough_beg)); %% ATTN: signal pust be all POSITIVE
                            temp1 = min(abs(data(trough_beg:trough_max)-side));

                            params(1).ind_fracpeak(ii) = find(abs(data(trough_beg:trough_max)-side)==temp1(1),1) + trough_beg-1;
                            params(1).eventsw(ii) = floor((params(1).ind_fracpeak(ii) - WOId(1))/size(data,1));
                            params(1).ind_insw(ii) = params(1).ind_fracpeak(ii) - double(params(1).eventsw(ii))*size(data,1);
                            params(1).risetime(ii) = trough_max-trough_beg; %% time from peak before trough to trough
                            params(1).max2min(ii) = abs(data(trough_max)-data(trough_beg)); %% amplitude of cycle from maximum before osc to trough
                            params(1).max(ii) =  data(trough_beg);
                            %% waveforms
                            %                  events(ii,:) = data(params(1).ind_fracpeak(ii)-WOI(1):params(1).ind_fracpeak(ii)+WOI(2));
                            %                  eventtrigger(ii,:) = intradata(params(1).ind_fracpeak(ii)-WOI(1):params(1).ind_fracpeak(ii)+WOI(2));



                            if bdebug
                                %                         plot(data(max(1,ind_xthres(i)-WOId(1)*2):ind_xthres(i)+WOId(2)*2))
                                %                 plot(abs(data(ind_xthres(i)-trough_beg:ind_xthres(i)+trough_max)-side))
                                hold on
                                %                 plot((ddata(ind_xthres(i)-WOId:ind_xthres(i))-mean(mean(ddata,1)))*100)

%                                 line([trough_beg-(ind_xthres(i)-WOId(1)) trough_beg-(ind_xthres(i)-WOId(1))], [min(data(ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2))) max(data(ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2)))])
%                                 line([trough_max-(ind_xthres(i)-WOId(1)) trough_max-(ind_xthres(i)-WOId(1))], [min(data(ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2))) max(data(ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2)))],'Color','r')
%                                 line([params(1).ind_fracpeak(ii)-(ind_xthres(i)-WOId(1)) params(1).ind_fracpeak(ii)-(ind_xthres(i)-WOId(1))], [min(data(ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2))) max(data(ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2)))],'Color','g')
                                line([trough_beg-(ind_xthres(i)) trough_beg-(ind_xthres(i))], [min(data(ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2))) max(data(ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2)))])
                                line([trough_max-(ind_xthres(i)) trough_max-(ind_xthres(i))], [min(data(ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2))) max(data(ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2)))],'Color','r')
                                line([params(1).ind_fracpeak(ii)-(ind_xthres(i)) params(1).ind_fracpeak(ii)-(ind_xthres(i))], [min(data(ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2))) max(data(ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2)))],'Color','g')
                                title(num2str(params(1).max2min(ii),'%d'));
                                ii
                                                        pause;
                            end

                            ii=ii+1;
                        else
                            params(1).skipped(i) = 1;
                        end% trough_beg~= trough_max
                    else
                        params(1).skipped(i) = 2;
                    end
                else
                    params(1).skipped(i) = 3;
                end
            else
                params(1).skipped(i) = 4;
            end
        else
            params(1).skipped(i) = 5;
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
s = sum(params(1).skipped>0);
params(1).ind_insw = params(1).ind_insw(1:end-s); %% MEM make into int32
params(1).risetime = params(1).risetime(1:end-s) ;
params(1).eventsw = params(1).eventsw(1:end-s) ; %% MEM make into int16
params(1).max2min = params(1).max2min(1:end-s);
params(1).ind_fracpeak = params(1).ind_fracpeak(1:end-s);
params(1).max = params(1).max(1:end-s);

