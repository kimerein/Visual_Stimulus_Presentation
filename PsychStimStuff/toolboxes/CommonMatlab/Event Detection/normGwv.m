function data = normGwv(data)
% function data = normGwv(data)
% normalize data so that trough is at -1
data = detrend(data);
data = -1*data./min(data);