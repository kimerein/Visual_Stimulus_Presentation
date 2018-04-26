function selev_ind = rmevents(ev_ind,minind)
% function selev_ind = rmevents(ev_ind,minind)
% remove events that in less than minind from previous event
% 
% note: doesn't a/c for case where previous event is removed and leaving a
% longer interval between event and pre-previous event
%
% BA030206

selev_ind = ev_ind([1 (find(diff(ev_ind) > minind) +1)']);
