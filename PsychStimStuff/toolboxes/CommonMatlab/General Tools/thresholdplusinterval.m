function sel_interval = thresholdplusinterval(data,thresh,min_ind)
% function sel_interval = thresholdplusinterval(data,thresh,min_ind)
% for each ROW of data (trial/sweep) returns intervals where data exceeds a value thresh.
% if min_ind (optional) specified, will return only intervals that are greater than min_ind
%  min_ind should be in units of indices of data
% USAGE: used to find peroids when powerspectrum is greather than threshold
% for LFP AMP vs IEI analysis
%
% BA053108

sel_interval  = [];
for i = 1:size(data,1)
    try

        th_data = data(i,:)>thresh; % find data
        start_ind = find( diff(th_data)==1);
        stop_ind = find( diff(th_data)==-1);
        if ~isempty(stop_ind) && ~isempty(start_ind)
            if stop_ind(1)<start_ind(1) % case of no first start
                start_ind = [1 start_ind];
                %         sel_interval = [sel_interval; i 1 stop_ind(1)];
                %         stop_ind = stop_ind(2:end);
            end

            if start_ind(end) > stop_ind(end)% case of no last stop
                %         flag= [i stop_ind(end) size(data,2)];
                stop_ind = [stop_ind size(data,2)];
            end

            for j=1:length(start_ind)
                sel_interval = [sel_interval; i start_ind(j) stop_ind(j)];
            end
        end
    catch
        keyboard
    end

end

if nargin>2 & ~isempty(sel_interval)% exclude below min_ind
    ind = find( sel_interval(:,3) - sel_interval(:,2) > min_ind);
    sel_interval = sel_interval(ind,:,:);
end

