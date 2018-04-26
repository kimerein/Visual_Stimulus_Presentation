%% USER input event detection


%% PARAMTERS
thres = -20; % in pA above/below baseline
sfilename = 'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\2006_05_08_0003sub.abf';
%%load file
output = importStruct_abf(sfilename,-1,0);
I = output.data(:,2:20);
swlength = length(I(:,1));
ii=1; jj=1; clear reject; accept = [];
for sw=1:size(I,2)
    If(:,sw) = filtfilt(B,A,double(I(:,sw)));
    %% remove mode of trace (assume mostly at baseline)
    [a x] = hist(If(:,sw),100);
    temp = max(a);
    temp = x(find(a == temp(1)));

    If(:,sw) = If(:,sw)-temp;

    I(:,sw) = I(:,sw)-temp;

    %     if bplot
    %         figure(1)
    %         clf
    %         plot(I(:,i),'-r');
    %         hold on;
    %         plot(If(:,i),'-b');
    %     end

    ind_bxthres = find(diff(If(:,sw)<thres)==1)+1;
    for i=1:length(ind_bxthres)

        figure(2) %% plot window around event 
        clf
        b =max([(ind_bxthres(i) - 50/1000*Ts) 1]);
        %         b =1;
        e = min ([(ind_bxthres(i) + 100/1000*Ts) size(I,1)]);
        plot(I(b:e,sw),'-r');
        hold on
        plot(If(b:e,sw),'-b','Linewidth',2);
        plot(ind_bxthres(i)-b,10,'.k','MarkerSize', 15);
        axis tight;
        title(['Sw: ' num2str(sw,'%d') 'Events: ' num2str(size(accept,2),'%d')])
        if jj>1
            if accept(jj-1,2) > b & accept(jj-1,1) == sw
                plot(accept(jj-1,2)-b,10,'.b','MarkerSize', 15);
            end
        end

        %% USER DECIDES to accept/reject
        reply = [];
        while isempty(reply)
            reply = input('Accept(a) Reject(r)  Quit(q):','s');
        end
        % reply = 'a';
        switch(reply)
            case 'r'
                reject(ii,:) = [sw ind_bxthres(i) -1];
                ii = ii+1;
            case 'a'
                %% find Xpercent of peak
                accept(jj,:) = [sw ind_bxthres(i)]
                jj = jj +1;
            case 'q'
                break;
        end


    end
end

da = accept;
dr = reject;

%% get selected events
WOI = [5 50]/1000*Ts; % window to take around event
lengtht = floor(WOI(2)-WOI(1));
evW = []; evWf = [];
for i = 1:size(da,1)
    b = da(i,2)-WOI(1);
    e = da(i,2)+WOI(2);
   if  (b>0)& (e< swlength)
       b = floor(b);
       e = b + lengtht;
       evW = [evW I(b:e,da(i,1))];%event waveform
       evWf = [evWf If(b:e,da(i,1))];%event waveform
   end
end
evW = evW';
evWf = evWf';

% baseline events
a =  repmat(mean(evW(:,1:WOI(1)/2),2),1,size(evW,2));
evW = evW - repmat(mean(evWf(:,1:WOI(1)*4/5),2),1,size(evW,2));
% mean
mevW = mean(evW);


%% PLOTTING
figure(3)
clf
plot([1:length(evW)]/Ts*1000,evW','-b')
hold on
plot([1:length(evW)]/Ts*1000,mevW,'-k','Linewidth',2)
axis tight
xlim([0 25])
xlabel('ms')
ylabel('pA')

   
%% Save
save('EPSPtimecourse', 'sfilename', 'thres', 'evW','evWf')

%% Fit
% Rise and Fall
% pick part of rising phase to fit
[rx y] = ginput(2)


y = mevW';
dt = 1/Ts;
riseWindow = round(rx/1000*Ts);


%% Fit by eye
tr = 2.5;
td =2.7;
%% Check fit
x = [0:.1:20]
figure(10)
clf
yf = exp(x/-td)-exp(x/-tr);
plot(x,yf/max(yf),'-k')
hold on
yy = y(riseWindow(1):riseWindow(2));
yy = (yy/max(yy))
yy = yy - min(yy)
yy = yy/max(yy)
xx = [1:size(yy,1)]./Ts.*1000;
plot(xx,yy)

