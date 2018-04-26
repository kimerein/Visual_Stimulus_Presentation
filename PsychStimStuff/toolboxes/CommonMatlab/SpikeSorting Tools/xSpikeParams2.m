function param = xSpikeParams2(data,dt,bsmooth)
% function param = xSpikeParams2(data,dt,bsmooth)
%% should be more robust then xSpikeParams
% dt = time between samples 
% Spike Parameters
% amplitude of P1,T1,P2
% width of P1,T1,P2
%% assume DC component is removed.
%% assume spikes take form Peak(P1),Trough(T1),Peak(P2).
%% parm = [P1 T1 P2 P1W T1W P2W T1F T1R P2R P2F T1Time ind_t1 ind_p2];
data = single(data);
data = data - mean(data,2)*ones(1,size(data,2));

if nargin<2
    dt =10.*ones(1,size(data,1));
end
if nargin <3
    bsmooth =1;
end
debug = 0;
%% Find Peak1
if bsmooth
    for i =1:size(data,1)
        data(i,:) = smooth(data(i,:),round(size(data,2)/13),'rloess');
        %     data(i,:) = smooth(data(i,:),'sgolay');
        if(mod(i,100) == 0)
            i
        end
    end
end
%% oversample to get accurate values for parameters
x = [1:size(data,2)];
step = .05; xs = [1:step:size(data,2)];
% tic
data = pchip(x,data,xs); %% if y is array size(y,2) == size(x)
% toc
% if debug
%     plot(x,data,'-b',xs,y,'-g');
% end

%% because freaking find doesnot find a [A B C D] in columns [1 2 3 4]
%% respectively.
param = zeros(size(data,1),12); 
for i = 1:size(data,1)
      try

%     ind_t1 = find(data(i,.125*size(data,2):end) == min(data(i,.125*size(data,2):end),1,'first')) +.125*size(data,2)-1; %% inserted offset because sometimes there is a transient?? at the beginning of the spikewaveform
    ind_t1 = find(data(i,floor(.125*size(data,2)):end) == min(data(i,floor(.125*size(data,2)):end))) +floor(.125*size(data,2))-1; %% inserted offset because sometimes there is a transient?? at the beginning of the spikewaveform
    ind_p1= find(data(i,:) == max(data(i,1:ind_t1)),1,'first');
    % use derivative because there are sometime multiple peaks after trough
    derdata = diff(data);
    ind_p2= find(derdata(i,floor(ind_t1+.125*size(data,2)):end)<0,1,'first') + (floor(ind_t1+.125*size(data,2))-1) +1;
    if isempty(ind_p2)
        ind_p2 = size(data,2);
    end
    %amplitude of peaks
    P1 = data(i,ind_p1);
    P2 = data(i,ind_p2);
    T1 = data(i,ind_t1);

    %% zerocrossing based
    % ind_P1T1 = (ind_p1 -1 + find(data == min(abs(data(ind_p1:ind_t1))))) ;
    % P1W = ind_P1T1 -find(data == min(abs(data(1:ind_p1)))) ; % width of P1  %% might be better to do some sort of quadratic fit
    % ind_T1P2 = (ind_t1 -1 + find(data == min(abs(data(ind_t1:ind_p2))))) ;
    % T1W = ind_T1P1 - ind_T1P2;%width of T1
    % ind_P2endz = ind_p2 -1 find(data == min(abs(data(ind_p2:end)))) ;
    % P2W = ind_P2end - ind_T1P2;%width of P2
    %
    %
    % T1F = ind_t1 - ind_P1T1;%trough-fall
    % T1R = ind_T1P2 - ind_t1; %trough-rise
    % R2R = ind_p2 - ind_T1P2;
    % T2F = ind_P2end - ind_p2;



%     %% 50% of peak based on left side of P1 and right side of P2
%     ind_begP1 = find(abs(data(i,1:ind_p1)-P1*.5) == min(abs(data(i,1:ind_p1)-P1*.5)),1,'last');
%     ind_P1T1 = ind_p1 -1 + find(abs(data(i,ind_p1:ind_t1)) == min(abs(data(i,ind_p1:ind_t1))),1,'last') ;
%     P1W = ind_P1T1 - ind_begP1 ; % width of P1  %% might be better to do some sort of quadratic fit
%     ind_T1P2 = ind_t1 -1 + find(abs(data(i,ind_t1:ind_p2)) == min(abs(data(i,ind_t1:ind_p2))),1,'last') ;
%     T1W = ind_T1P2 -ind_P1T1;%width of T1
%     ind_P2end = ind_p2 -1 + find(abs(data(i,ind_p2:end)-P2*.5) == min(abs(data(i,ind_p2:end)-P2*.5)),1,'first') ;
%     P2W = ind_P2end - ind_T1P2;%width of P2

%     T1F = ind_t1 - ind_P1T1;%trough-fall
%     T1R = ind_T1P2 - ind_t1; %trough-rise
%     P2R = ind_p2 - ind_T1P2;
% %     P2F = ind_P2end - ind_p2;
% 
%    if debug  %% show parameters on plot
%         plot(data(i,:),'-g');
%         line([ind_t1 ind_t1],[min(data(i,:)) max(data(i,:))],'Color','r','LineWidth',0.5); %VERTICAL
%         line([ind_p1 ind_p1],[min(data(i,:)) max(data(i,:))],'Color','b','LineWidth',0.5); %VERTICAL
%         line([ind_p2 ind_p2],[min(data(i,:)) max(data(i,:))],'Color','k','LineWidth',0.5); %VERTICAL
% 
%         line([ind_begP1 ind_begP1],[min(data(i,:)) max(data(i,:))],'Color','g','LineWidth',0.5); %VERTICAL
%         line([ind_P1T1 ind_P1T1],[min(data(i,:)) max(data(i,:))],'Color','c','LineWidth',0.5); %VERTICAL
%         line([ind_T1P2 ind_T1P2],[min(data(i,:)) max(data(i,:))],'Color','y','LineWidth',0.5); %VERTICAL
%         line([ind_P2end ind_P2end],[min(data(i,:)) max(data(i,:))],'Color','m','LineWidth',0.5); %VERTICAL
%         pause
%     end

  %% FWHM for all
    ind_begP1 = find(abs(data(i,1:ind_p1)-P1*.5) == min(abs(data(i,1:ind_p1)-P1*.5)),1,'last');
    ind_P1T1 = ind_p1 -1 + find(abs(data(i,ind_p1:ind_t1)) == min(abs(data(i,ind_p1:ind_t1))),1,'last') ;
    % index of falling phase of P1
    ind_P1F = ind_p1 -1 + find(abs(data(i,ind_p1:ind_P1T1)-P1*.5) == min(abs(data(i,ind_p1:ind_P1T1)-P1*.5)),1,'last') ;
    P1W = ind_P1F - ind_begP1 ; % width of P1  %% might be better to do some sort of quadratic fit
    % index of falling phase of T1  AND rising
%     ind_T1F = ind_P1T1 -1 + find(abs(data(i,ind_P1T1:ind_t1)-T1*.5) == min(abs(data(i,ind_P1T1:ind_t1)-T1*.5)),1,'last') ;
%     ind_T1R = ind_t1 -1 + find(abs(data(i,ind_t1:ind_p2)-T1*.5) == min(abs(data(i,ind_t1:ind_p2)-T1*.5)),1,'last') ;
    ind_T1R = ind_t1 -1 + find((data(i,ind_t1:ind_p2)-T1*.5) > 0,1,'first') ;
    ind_T1F = ind_P1T1 -1 + find((data(i,ind_P1T1:ind_t1)-T1*.5) > 0,1,'last') ;

    T1W = ind_T1R -ind_T1F;%width of T1
      
    ind_T1P2 = ind_t1 -1 + find((data(i,ind_t1:ind_p2)) < 0,1,'last') ;
    
%     ind_P2R = ind_T1P2 -1 + find(abs(data(i,ind_T1P2:ind_p2)-P2*.5) == min(abs(data(i,ind_T1P2:ind_p2)-P2*.5)),1,'first') ;
   ind_P2R = ind_T1P2 -1 + find((data(i,ind_T1P2:ind_p2)-P2*.5) >0,1,'first') ;
    ind_P2F = ind_p2 -1 + find((data(i,ind_p2:end)-P2*.5) <0,1,'first') ;
    if isempty(ind_P2F)
        ind_P2F=size(data,2);
    end
    P2W = ind_P2F - ind_P2R;%width of P2

    T1F = ind_t1 - ind_T1F;%trough-fall
    T1R = ind_T1R - ind_t1; %trough-rise
    P2R = ind_p2 - ind_P2R;
    P2F = ind_P2F - ind_p2;

    if debug  %% show parameters on plot
        figure(1000)
        plot(data(i,:),'-g');
        line([ind_t1 ind_t1],[min(data(i,:)) max(data(i,:))],'Color','r','LineWidth',0.5); %VERTICAL
        line([ind_p1 ind_p1],[min(data(i,:)) max(data(i,:))],'Color','b','LineWidth',0.5); %VERTICAL
        line([ind_p2 ind_p2],[min(data(i,:)) max(data(i,:))],'Color','k','LineWidth',0.5); %VERTICAL

        line([ind_P2R ind_P2R],[min(data(i,:)) max(data(i,:))],'Color','g','LineWidth',0.5); %VERTICAL
        line([ind_P2F ind_P2F],[min(data(i,:)) max(data(i,:))],'Color','c','LineWidth',0.5); %VERTICAL
        line([ind_T1F ind_T1F],[min(data(i,:)) max(data(i,:))],'Color','y','LineWidth',0.5); %VERTICAL
        line([ind_T1R ind_T1R],[min(data(i,:)) max(data(i,:))],'Color','m','LineWidth',0.5); %VERTICAL
        pause
    end
        if size(dt,1) ==1
            param(i,:) = [P1 T1 P2 P1W*step*dt T1W*step*dt P2W*step*dt T1F*step*dt T1R*step*dt P2R*step*dt P2F*step*dt ind_t1*step*dt ind_p2*step*dt];
        else
            param(i,:) = [P1 T1 P2 P1W*step*dt(i,1) T1W*step*dt(i,1) P2W*step*dt(i,1) T1F*step*dt(i,1) T1R*step*dt(i,1) P2R*step*dt(i,1) P2F*step*dt(i,1) ind_t1*step*dt(i,1) ind_p2*step*dt(i,1)];
        end
  catch
      i
  end
end
