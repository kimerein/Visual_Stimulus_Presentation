function [ind_pkCont ind_trCont ind_pk ind_tr] = detPeakTrough(data,th)
% function [ind_pk ind tr] = detPeakTrough(data,th)
%
% returns indices of peaks and troughs in data,
% INPUT:
%        data  (data in vector)
%        th (optional) threshold for difference between peak and trough

bdebug = 0;
ind_pk = findpeaks(data);
ind_tr = findpeaks(-data);


if bdebug
    a=0;
    for i= 1:length(ind_pk)
        a = a +length(ind_pk(i).loc);
    end
    a
    figure(1);
    clf;
    i =1;
    plot(ind_pk(i).loc,0*ones(size(ind_pk(i).loc)),'.r')
    hold on
    plot(ind_pk(i).zerodiff,0*ones(size(ind_pk(i).zerodiff)),'.k')
    hold on
    plot(ind_tr(i).loc,0*ones(size(ind_tr(i).loc))/2,'.g')
    hold on
    plot(data(:,i)/max(data(:,i)));
end

for i= 1:length(ind_pk)
    %first peak must be before first trough.
    if ~isempty(ind_tr(i).loc) && ~isempty(ind_pk(i).loc)
        if ind_tr(i).loc(1)<ind_pk(i).loc(1);
            ind_tr(i).loc = ind_tr(i).loc(2:end);
        end
        %last peak must be before last trough.
        if ind_tr(i).loc(end)<ind_pk(i).loc(end);
            ind_pk(i).loc = ind_pk(i).loc(1:end-1);
        end

        if length(ind_pk(i).loc) ~= length(ind_tr(i).loc)
            error ('ind_tr and ind_pk are not the same length');
        end
    end
    if nargin >1 && ~isempty(th)
        tmp = find(abs(data(ind_tr(i).loc,i)-data(ind_pk(i).loc,i)) > abs(th));
        ind_tr(i).loc = ind_tr(i).loc(tmp);
        ind_pk(i).loc = ind_pk(i).loc(tmp);
    end
    if bdebug
        figure(6);clf;
        plot(data(:,i))
        hold on;plot(ind_pk(i).loc,zeros(size(ind_pk(i).loc)),'.r')
        hold on;plot(ind_tr(i).loc,-50*ones(size(ind_tr(i).loc)),'.g')
        plot(ind_pk(i).loc(tmp),zeros(size(ind_pk(i).loc(tmp))),'.k')
        plot(ind_tr(i).loc(tmp),-50*ones(size(ind_tr(i).loc(tmp))),'.k')
        plot(ind_tr(i).loc,abs(data(ind_tr(i).loc)-data(ind_pk(i).loc)),'k')
    end
end

% convert into one long vector for compatibility with how thresdata.m works
tmppk = [];tmptr = [];
for i= 1:length(ind_pk)
    tmppk = [tmppk; ind_pk(i).loc+(i-1)*size(data,1)];
    tmptr= [tmptr; ind_tr(i).loc+(i-1)*size(data,1)];
end

ind_pkCont = tmppk;
ind_trCont= tmptr;

if bdebug
    figure(1);
    clf;
    i =1;
    plot(ind_pk,0*ones(size(ind_pk)),'.r')
    hold on
    plot(ind_tr,0*ones(size(ind_tr)),'.g')
    hold on
    plot(data/max(data));
end
