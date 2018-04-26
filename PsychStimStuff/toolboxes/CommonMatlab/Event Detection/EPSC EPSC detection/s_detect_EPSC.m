%% EPSC detection
%% GOOD START 

DOESNOT QUITE WORK YET


%% PARAMTERS
thres = -5; % in pA above/below baseline
diff_thres = 0.005 % threshold below which the derivative of I is considered to be part of a peak
pkWin = 1 ; %(ms) a Window around peak for which derivative is allowed to be less then diff_thres
%% outside of this window algorithm decides that there is a 2nd peak so the
BEFORE_bxthreshold = 1; %%(ms) this paramter depends on the rise time for assume fast rise time

%%load file
output = importStruct_abf('C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\2006_05_08_0003sub.abf',-1,0);

I = output.data(:,2:20);
bplot = 0;
Ts = 1/output.dt;
CUT = 300;%Hz
[B,A] = butter(2,2*CUT/Ts,'low');

pkWin = pkWin/1000 * Ts; % convert into samples
BEFORE_bxthreshold = int32(BEFORE_bxthreshold/1000*Ts);
for i=1:size(I,2)
    If(:,i) = filtfilt(B,A,double(I(:,i)));
    %% remove mode of trace (assume mostly at baseline)
    [a x] = hist(If(:,i),100);
    temp = max(a);
    temp = x(find(a == temp(1)));

    If(:,i) = If(:,i)-temp;

    if bplot
        figure(1)
        clf
        plot(I(:,i),'-r');
        hold on;
        plot(If(:,i),'-b');
    end
end

% tempIf = (If<thres).*If)
if bplot
    figure(2)
    clf
    plot(I(:,2),'-r');
    hold on;
    hold on
    temp = find(ind_xthres < size(I,1)*2& ind_xthres > size(I,1))
    plot(ind_xthres(temp),zeros(1,size(temp,2)),'xb','MarkerSize', 10);
end

%% Check fo r
jj = 1; ii =1; clear accept; clear reject;
for sw = 1: size(I ,2)
    ind_bxthres = find(diff(If(:,sw)<thres)==1)+1;
    ind_exthres = find(diff(If(:,sw)<thres)==-1)+1;
    % check that bxthes and exthres match
    if  ind_bxthres(1) > ind_exthres(1)
        ind_exthres = ind_exthres(2:end)
    end
    if  ind_bxthres(end) > ind_exthres(end)
        ind_bxthres = ind_bxthres(1:end-1)
    end
    ind_bxthres = int32(ind_bxthres);
    ind_exthres = int32(ind_exthres);
% %% DEBUG
% plot(If(:,sw))
% hold on
% plot(ind_bxthres,zeros(1,length(ind_bxthres)),'.r')
% plot(ind_exthres,zeros(1,length(ind_exthres)),'.b')
%%
    for i=1:length(ind_bxthres)
        %% THRESHOLD EVENTS
        eventAmp = max(If(ind_bxthres(i):ind_exthres(i),sw));
        ind_peak = int32(find(If(ind_bxthres(i):ind_exthres(i),sw) == eventAmp)) +ind_bxthres(i) -1 ;

        % time when derivative is close to 0 in indexes
        temp = int32(find(abs(diff(If(ind_bxthres(i):ind_exthres(i),sw)))< diff_thres)) +ind_bxthres(i) -1 ; %% time of peaks from derivative

        % find derivative is low just before the peak (used for baseline adjusment
        if (ind_bxthres(i)- BEFORE_bxthreshold) > 0
            temp2 = int32(find(abs(diff(If((ind_bxthres(i)- BEFORE_bxthreshold):ind_bxthres(i),sw)))< diff_thres,1,'last')) +ind_bxthres(i)- BEFORE_bxthreshold ;
            if isempty(temp2)
                temp2 = ind_bxthres(i)-1/1000*Ts;
            end
            
            if pkWin<max(ind_peak*ones(length(temp),1,'int32')-temp) %% check there deriviative is close to 0 for NO MORE than 1ms around the peak
                reject(ii,:) = [sw ind_bxthres(i) eventAmp -1];
                ii = ii+1;
            else
                %% find Xpercent of peak
                accept(jj,:) = [sw temp2 ind_bxthres(i) eventAmp]
                jj = jj +1;
            end
        else
            reject(ii,:) = [sw ind_bxthres(i) eventAmp -2];
        end
    end
end
figure(5)
clf
i=1
plot(If(:,i));
hold on
ind = find (accept(:,1)==i);
plot(accept(ind,2),zeros(1,length(ind)),'.r')
ind = find (reject(:,1)==i);
plot(reject(ind,2),zeros(1,length(ind)),'.k')

%
% temp = diff((If<thres).*If);
% plot(temp(:,1),'-k')

%%