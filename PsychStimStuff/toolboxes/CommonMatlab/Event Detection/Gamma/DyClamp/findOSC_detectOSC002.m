%% Detect oscillations in prodata(lfpdata filtered)
%% extract intradata triggered by detected oscillations
%% changed to use 20% peak as aligning point rather then threshold crossing
bdebug = 0;
%% tempoarary
prodata =   filtfilt(B,A,mondata);
xthres = find(diff(prodata < thres)==1);
dprodata  = diff(prodata);

%% pre process prodata
temp = min(min(prodata));
if temp < 0 %% make all postive if not;
    prodata2 = prodata - temp*1.1;
end

%% Extract event
WOI = round((1/output.dt/1000)).*[40 40];
WOId = round((1/output.dt/1000)).*20% (ms) too look for first positive derivitive
pkWOI = round((1/output.dt/1000).*[4 6]); %% for peak search (must be small then WOI)

events = zeros(size(xthres,1),sum(WOI)+1); eventtrigger = zeros(size(xthres,1),sum(WOI)+1);
eventsw = zeros(1,size(xthres,1));risetime = eventsw;
timeinsw = eventsw;
max2min = eventsw;

Ncycles = size(xthres,1);
ii = 1;
for i = 1: Ncycles
    try

        %         check same sweep
        if floor((xthres(i) - WOI(1))/size(prodata2,1)) == floor((xthres(i) + WOI(2))/size(prodata2,1))
            %% find where the trough begins by finding were the derivative goes positive
            trough_beg = find(dprodata(xthres(i)-WOId:xthres(i)) >0,1,'last') + xthres(i)-WOId; %% add 1 because of diff
            %% find trough
            if ~isempty(trough_beg) %% no rise trough is found in dWOI.. rise time is too slow
                trough_max = find(dprodata(xthres(i):xthres(i)+round(WOId)) > 0,1,'first') + xthres(i); %% add 1 because of diff
             if ~isempty(trough_max) %% no trough is found in dWOI.. rise time is too slow
                 %% find amp of 20% of peak
                 side = prodata2(trough_beg)-0.2*abs(prodata2(trough_max)-prodata2(trough_beg)); %% ATTN: signal pust be all POSITIVE
                 temp1 = min(abs(prodata2(trough_beg:trough_max)-side));
                 ind_sid = find(abs(prodata2(trough_beg:trough_max)-side)==temp1(1)) + trough_beg-1;
                 eventsw(ii) = floor((ind_sid - WOI(1))/size(prodata2,1));
                 timeinsw(ii) = ind_sid - double(eventsw(ii))*size(prodata2,1);
                 risetime(ii) = trough_max-trough_beg; %% time from peak before trough to trough
                 max2min(ii) = abs(prodata2(trough_max)-prodata2(trough_beg)); %% amplitude of cycle from maximum before osc to trough
                 
                 events(ii,:) = prodata2(ind_sid-WOI(1):ind_sid+WOI(2));
                 eventtrigger(ii,:) = intradata(ind_sid-WOI(1):ind_sid+WOI(2));
                 ii=ii+1;

                 if bdebug
                     figure(7)
                     clf
                     plot(prodata2(xthres(i)-WOI(1):xthres(i)+WOI(2)))
                     %                 plot(abs(prodata2(xthres(i)-trough_beg:xthres(i)+trough_max)-side))
                     hold on
                     %                 plot((dprodata(xthres(i)-WOId:xthres(i))-mean(mean(dprodata,1)))*100)

                     line([trough_beg-(xthres(i)-WOI(1)) trough_beg-(xthres(i)-WOI(1))], [min(prodata2(xthres(i)-WOI(1):xthres(i)+WOI(2))) max(prodata2(xthres(i)-WOI(1):xthres(i)+WOI(2)))])
                     line([trough_max-(xthres(i)-WOI(1)) trough_max-(xthres(i)-WOI(1))], [min(prodata2(xthres(i)-WOI(1):xthres(i)+WOI(2))) max(prodata2(xthres(i)-WOI(1):xthres(i)+WOI(2)))],'Color','r')
                     line([ind_sid-(xthres(i)-WOI(1)) ind_sid-(xthres(i)-WOI(1))], [min(prodata2(xthres(i)-WOI(1):xthres(i)+WOI(2))) max(prodata2(xthres(i)-WOI(1):xthres(i)+WOI(2)))],'Color','g')
                     title(num2str(risetime(ii)*output.dt*1000,'%1.1f'));
                     ii
                     %                pause;
                 end
             end
            end
        end

        %         plot(events(i,:))
        %         pause;
        if mod(i,50) ==0 %% Print progress
            sprintf([num2str(i) ' of ' num2str(Ncycles)]);
        end
    catch
        sprintf(['ERROR at i=' num2str(i,'%d') ' of ' num2str(Ncycles)])
    end
end

% %% find ISI and apply range filter
% ISI = diff(xthres)*output.dt*1000;
% Strough= xthres(find ( (ISI>=ISIrange(2)).*(ISI<=ISIrange(1)))+1);
% % Strough = trough;
% clear trough;
% 
% eventsw = int16(eventsw); %% to save space