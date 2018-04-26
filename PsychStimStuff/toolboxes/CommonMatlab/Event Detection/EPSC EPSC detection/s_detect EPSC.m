%% EPSC detection

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
BEFORE_bxthreshold = BEFORE_bxthreshold/1000*Ts;
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

%% THRESHOLD EVENTS
ind_bxthres = find(diff(If<thres)==1)+1;
ind_exthres = find(diff(If<thres)==-1)+1;
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
   jj = 1; ii =1;
   for i=1:size(ind_bxthres)
       
       eventAmp = max(If(ind_bxthres(i):ind_exthres(i)));
       ind_peak = find(If(ind_bxthres(i):ind_exthres(i)) == eventAmp(i)) +ind_bxthres(i) -1 ;
       
       % time when derivative is close to 0 in indexes
       temp = find(diff(If(ind_bxthres(i):ind_exthres(i)))< diff_thres) +ind_bxthres(i) -1 ; %% time of peaks from derivative

       % find derivative is low just before the peak (used for baseline adjusment
       temp2 = find(diff(If((ind_bxthres(i)- BEFORE_bxthreshold):ind_bxthres(i))< diff_thres,1,'last') +ind_bxthres(i)- BEFORE_bxthreshold ;       

       if pkWin<max(ind_peak*ones(length(temp),1)-temp) %% check there deriviative is close to 0 for NO MORE than 1ms around the peak
           reject(ii) = i;
           ii = ii+1;
       else
           %% find Xpercent of peak
           accept(jj,:) = [temp2 ind_bxthres(i) eventAmp]
           jj = jj +1;
       end
   end
      
    figure(5)
    clf
plot(If(1:end));
hold on
plot(accept(1,:),zeros(1,length(accept)),'.r')

% 
% temp = diff((If<thres).*If);
% plot(temp(:,1),'-k')

%% 