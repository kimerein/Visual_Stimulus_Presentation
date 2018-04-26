%% Detect oscillations in prodata(lfpdata filtered)
%% extract intradata triggered by detected oscillations
%% changed to use 20% peak as aligning point rather then threshold crossing
a= prodata < thres;
a = [diff(prodata); zeros(1,size(prodata,2))].*a;
j = 1;

%% looks like it is aligned to threshold crossing 
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
bdebug = 0;
%% Extract event
WOI = round((1/output.dt/1000)).*[40 10];
pkWOI = round((1/output.dt/1000).*[4 6]); %% for peak search (must be small then WOI)
events = zeros(size(Strough,1),sum(WOI)+1); eventtrigger = zeros(size(Strough,1),sum(WOI)+1);
Ncycles = size(Strough,2);
ii = 1;

%% pre process prodata
temp = min(min(prodata));
if temp < 0 %% make all postive if not;
    prodata = prodata - temp*1.1;
end
for i = 1: Ncycles
    try

        %         check same sweep
        if floor((Strough(i) - WOI(1))/size(prodata,1)) == floor((Strough(i) + WOI(2))/size(prodata,1))
            %% find 20% of peak
            side = (max(prodata(Strough(i)-pkWOI(1):Strough(i)+pkWOI(2)))-min(prodata(Strough(i)-pkWOI(1):Strough(i)+pkWOI(2))))*.20;
            temp1 = min(abs(prodata(Strough(i)-pkWOI(1):Strough(i)+pkWOI(2))-side));
            ind_sid = find(abs(prodata(Strough(i)-pkWOI(1):Strough(i)+pkWOI(2))-side)==temp1(1)) + Strough(i)-pkWOI(1)-1;
            eventsw(ii) = floor((ind_sid - WOI(1))/size(prodata,1));
            timeinsw(ii) = ind_sid - (eventsw(ii))*size(prodata,1);

            events(ii,:) = prodata(ind_sid-WOI(1):ind_sid+WOI(2));
            eventtrigger(ii,:) = intradata(ind_sid-WOI(1):ind_sid+WOI(2));
            ii=ii+1;
            
            if bdebug
               figure(7)
               plot(prodata(ind_sid-WOI(1):ind_sid+WOI(2)))
               ii
%                pause;
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
eventsw = int16(eventsw); %% to save space