%% Detect oscillations in prodata(lfpdata filtered)
%% extract intradata triggered by detected oscillations
a= prodata < thres;
a = [diff(prodata); zeros(1,size(prodata,2))].*a;
j = 1;
    for (i = 2: size(a,1)*size(a,2))
        if a(i) > 0 & a(i-1) <= 0 % if trough
            if min(abs(a(i)),abs(a(i-1)))==abs(a(i))
                trough(j) = i;
                baug = 1;
            else
                trough(j) = i-1;
                baug = 1;
            end
            if baug
                j = j+1;
                %% % plot to debug detection of troughs
                debug = 0;
                if debug
                    figure(100);
                    plot(a(1:19999),'k') %%% DATA SPECIFIC
                    hold on
                    plot(prodata(1:19999),'b')
                    plot(trough(j),.5,'.r')
                    pause;
                end
            end
        end
    end
clear a

%% find ISI and apply range filter
ISI = diff(trough)*output.dt*1000;
Strough= trough(find ( (ISI>=ISIrange(2)).*(ISI<=ISIrange(1)))+1);
% Strough = trough;
clear trough;

%% Extract event
WOI = (1/output.dt/1000).*[40 10];
events = zeros(size(Strough,1),sum(WOI)+1); eventtrigger = zeros(size(Strough,1),sum(WOI)+1); 
Ncycles = size(Strough,2);
for i = 1: Ncycles
    try
        %         check same sweep
        if int32((Strough(i) - WOI(1))/size(prodata,1)) == int32((Strough(i) + WOI(2))/size(prodata,1))
            eventsw = floor((Strough(i) - WOI(1))/size(prodata,1));
            timeinsw(i) = Strough(i) - (eventsw)*size(prodata,1);
            events(i,:) = prodata(Strough(i)-WOI(1):Strough(i)+WOI(2));
            eventtrigger(i,:) = intradata(Strough(i)-WOI(1):Strough(i)+WOI(2));
        end
        %         plot(events(i,:))
        %         pause;

    catch
        i
    end
end