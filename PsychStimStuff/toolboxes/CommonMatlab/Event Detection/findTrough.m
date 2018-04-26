function [params] = findTrough(data,datatype,ind_xthres,WOId)
% function [params] = findTrough(data,datatype,ind_xthres,WOId)
%
% Finds index where amplitude is minimum (trough), searching in the window specified
% by WOId
%
% INPUT: data  = data should contain 1 channel data from an experiment
%                each row should contain 1 sweep
%      NOTE: data should be filtered because diff will be used to find
%      peaks and troughs in data
%        datatype = 1 if gapfree, 0 if episodic
%
%       ind_xthres = index of threshold crossing (from thresh function
%       WOId = 1x2 array specifying the window to look for a trough
%                     [before after]
% OUTPUT:
%         params(1).ind_tr = index of trough
%             if WOId does not exist for a ind_xthres (index in experiment)
%         params  = struct containing
%               params(1).ind_insw =  index in sweep
%               params(1).eventsw = sweep number; %% MEM make into int16
%               params(1).trAmp = trough amplitude
%               params(1).skipped  = 1 if is skipped
% BA011607
bdebug =0;


Ncycles = length(ind_xthres);


%% DECLARE PARAMS
params = struct([])
params(1).ind_insw = zeros(1,Ncycles); %% MEM make into int32
params(1).eventsw = zeros(1,Ncycles); %% MEM make into int16
params(1).ind_tr = params(1).ind_insw;
params(1).trAmp = NaN(1,Ncycles);
params(1).skipped = zeros(1,Ncycles,'int16')
ii = 1;
for i = 1: Ncycles
    try
        %         check same sweep
        samesw = 0;
        if datatype
            samesw = 1;
        elseif floor((ind_xthres(i) - WOId(1))/size(data,1)) == floor((ind_xthres(i) + WOId(2))/size(data,1))
            samesw = 1;
        end
        if samesw
            %% find trough
            params(1).trAmp(ii) = min(data(ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2)));
            params(1).ind_tr(ii) = find(data(ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2))==params(1).trAmp(ii),1,'first') + ind_xthres(i)-WOId(1)-1; %% add 1 because of diff
            params(1).eventsw(ii) = floor((params(1).ind_tr(ii) - WOId(1))/size(data,1));
            params(1).ind_insw(ii) = params(1).ind_tr(ii) - double(params(1).eventsw(ii))*size(data,1);
            %% waveforms
            %             events(ii,:) = data(params(1).ind_tr(ii)-WOI(1):params(1).ind_tr(ii)+WOI(2));
            %% DEBUG
            if bdebug
                figure(7)
                clf
                plot(data(ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2)))
                %                         plot(data(max(1,ind_xthres(i)-WOId(1)*2):ind_xthres(i)+WOId(2)*2))
                %                 plot(abs(data(ind_xthres(i)-trough_beg:ind_xthres(i)+trough_max)-side))
                hold on
                %                 plot((ddata(ind_xthres(i)-WOId:ind_xthres(i))-mean(mean(ddata,1)))*100)

                line([params(1).ind_tr(ii)-(ind_xthres(i)-WOId(1)) params(1).ind_tr(ii)-(ind_xthres(i)-WOId(1))], [min(data(ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2))) max(data(ind_xthres(i)-WOId(1):ind_xthres(i)+WOId(2)))],'Color','g')
                title(num2str(params(1).trAmp(ii),'%d'));
                ii
                %                         pause;
            end

            ii=ii+1;
        else
            params(1).skipped(i) = 1;
        end% trough_beg~= trough_max



        %% Print progress
        if mod(i,50) ==0
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
params(1).eventsw = params(1).eventsw(1:end-s) ; %% MEM make into int16
params(1).ind_tr = params(1).ind_tr(1:end-s);
params(1).trAmp = params(1).trAmp(1:end-s);
