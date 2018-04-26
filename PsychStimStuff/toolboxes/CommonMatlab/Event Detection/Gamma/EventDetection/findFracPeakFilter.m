function [params] = findFracPeakFilter(data, dt, dfilter, ind_xthres,frac,WOId)
% function [params] = findFracPeakFilter(data, dt, dfilter, ind_xthres,frac,WOId)
%
% Similar to findFracPeak (see documentation for more details), but input 'data' can be unfiltered and is filtered according to dfilter.
% ind_xthres if found using filtered 'data' but amplitude is take from unfiltered.`
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
%       dfilter =[HIGH_PASS_Hz LOW_PASS_Hz  boolean_60Hz_removal]
%                   e.g. [20 200 0]
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
bdebug =0;

display('NOTE: data must be supplied such that each ROW is a new sweep/channel')

df = prepdata (data,0,dt,dfilter)'; %preprocess input each sweep is a ROW transposed so df is COL
if size(data,1)>1 % must do this for max2min to work on data with more than one sweep
    data = data';
    data = data(1:end);
end
% if isrowvector(df)
%     df = df';
% end
ddf = diff(df,1,1); % takes derivative along COL (each COL should be a sweep/episode)

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
for i = 1: Ncycles
    try
% ind_xthres(i)/length(df)
%         if floor((ind_xthres(i) - WOId(1))/size(df,1)) == floor((ind_xthres(i) + WOId(2))/size(df,1))  %  
           temp_sw = floor((ind_xthres(i) - WOId(1)-1)/size(df,1));
        if temp_sw == floor((ind_xthres(i) + WOId(2))/size(df,1))  %         check same sweep
            %% find where the trough begins by finding were the derivative goes positive
            if ind_xthres(i)> WOId(1) %check first point as WOId(1) before it
                trough_beg = find(ddf((ind_xthres(i)-WOId(1):ind_xthres(i))-(temp_sw)*size(df,1),temp_sw+1)>0.0,1,'last') + ind_xthres(i)-WOId(1)-1; %% add 1 because of diff
                %% find trough
                if bdebug
                    figure(7)
                    clf
%                     plot([-WOId(1):0],ddf(ind_xthres(i)-WOId(1):ind_xthres(i)) >0);  
                    hold on;
                    plot([-WOId(1):0],ddf((ind_xthres(i)-WOId(1):ind_xthres(i))-(temp_sw)*size(df,1),temp_sw+1))
                    plot([-WOId(1):WOId(2)],df((ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2))-(temp_sw)*size(df,1),temp_sw+1),'-k');

%                    plot([-WOId(1):WOId(2)],ddf((ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2))-(temp_sw)*size(df,1),temp_sw+1),'-k');

% %                     hold on
% %                   (BROKEN)  plot([-WOId(1):WOId(2)],data((ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2))-(temp_sw)*size(df,1),temp_sw+1));
%                              PLOT entire sweep               figure(7)
% %                     clf
% %                     plot(df(:,temp_sw+1),'-k');
% %                     hold on
% %                                         plot(ddf(:,temp_sw+1),'-k');
% %                                       plot(data(:,temp_sw+1),'-k');
% 
%                     xlim([ind_xthres(i)-WOId(1) ind_xthres(i)+WOId(2)]-(temp_sw)*size(df,1))
%                     ylim([min(data((ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2))-(temp_sw)*size(df,1),temp_sw+1)) max(data((ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2))-(temp_sw)*size(df,1),temp_sw+1))]);
                    
                    if ind_xthres(i)-WOId(1)>=1
                        temp = ind_xthres(i)-WOId(1)*10;
                        if temp >0
%                             plot([-WOId(1)*10:WOId(2)*10],df((temp:ind_xthres(i)+WOId(2)*10) -(temp_sw)*size(df,1),temp_sw+1));
                        end
                    end
                end
                if ~isempty(trough_beg) %% no rise trough is found in dWOI.. rise time is too slow
%                     trough_max = find(ddf(trough_beg:ind_xthres(i)+round(WOId(2))) > 0,1,'first') + trough_beg; %% add 1 because of diff
                                      trough_max = find(ddf((trough_beg+2:ind_xthres(i)+round(WOId(2)))-temp_sw*size(df,1),temp_sw+1)>0.0,1,'first') + trough_beg+1;
 
                    if ~isempty(trough_max) %% no trough is found in dWOI.. rise time is too slow
                        %% find amp of frac of peak
                        if trough_beg~= trough_max
                            side = df(trough_beg)-frac*abs(df(trough_max)-df(trough_beg)); %% ATTN: signal must be all POSITIVE
                            temp1 = min(abs(df(trough_beg:trough_max)-side));

                            params(1).ind_fracpeak(ii) = find(abs(df(trough_beg:trough_max)-side)==temp1(1),1) + trough_beg-1;
                            params(1).eventsw(ii) = floor((params(1).ind_fracpeak(ii) - WOId(1))/size(df,1))+1;
                            params(1).ind_insw(ii) = params(1).ind_fracpeak(ii) - double(params(1).eventsw(ii))*size(df,2);
                            params(1).risetime(ii) = trough_max-trough_beg; %% time from peak before trough to trough
                            dataWIN = [1 1]; %ms
                            params(1).max(ii) =  max(data((round((trough_beg+1)- dataWIN(1)/1000/dt)):round((trough_beg+1)+dataWIN(1)/1000/dt)));
                            params(1).max2min(ii) = range(data(round(((trough_beg+1)):round((trough_max+1)))));%% amplitude of cycle from maximum before osc to trough
                            %% waveforms
                            %                  events(ii,:) = df(params(1).ind_fracpeak(ii)-WOI(1):params(1).ind_fracpeak(ii)+WOI(2));
                            %                  eventtrigger(ii,:) = intradf(params(1).ind_fracpeak(ii)-WOI(1):params(1).ind_fracpeak(ii)+WOI(2));



                            if bdebug
                                %                         plot(df(max(1,ind_xthres(i)-WOId(1)*2):ind_xthres(i)+WOId(2)*2))
%                                                 plot(abs(df(max(ind_xthres(i)-trough_beg,1):ind_xthres(i)+trough_max)-side))
%                                 hold on
%                                                 plot((ddf(ind_xthres(i)-WOId:ind_xthres(i))-mean(mean(ddf,1)))*100)
% 
%                                 line([trough_beg-(ind_xthres(i)-WOId(1)) trough_beg-(ind_xthres(i)-WOId(1))], [min(df(ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2))) max(df(ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2)))])
%                                 line([trough_max-(ind_xthres(i)-WOId(1)) trough_max-(ind_xthres(i)-WOId(1))], [min(df(ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2))) max(df(ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2)))],'Color','r')
%                                 line([params(1).ind_fracpeak(ii)-(ind_xthres(i)-WOId(1)) params(1).ind_fracpeak(ii)-(ind_xthres(i)-WOId(1))], [min(df(ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2))) max(df(ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2)))],'Color','g')
                                line([trough_max-(ind_xthres(i)) trough_max-(ind_xthres(i))], [min(df(ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2))) max(df(ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2)))],'Color','r')
                                line([trough_beg-(ind_xthres(i)) trough_beg-(ind_xthres(i))], [min(df(ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2))) max(df(ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2)))])
                                line([params(1).ind_fracpeak(ii)-(ind_xthres(i)) params(1).ind_fracpeak(ii)-(ind_xthres(i))], [min(df(ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2))) max(df(ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2)))],'Color','g')
%                                  title(num2str(params(1).risetime(ii),'%d'));
                                title(num2str(params(1).max2min(ii),'%d'));
%                                 temp = (temp_sw*length(df));
%                                  line([trough_max trough_max], [min(df(ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2))) max(df(ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2)))],'Color','r')
                              ii;
%                                                         pause;
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

