function f = ErrRiseTimefit(D)
% function ErrRiseTimefit(D)
%
% error function for fitting PSC with correct 10-90 rise time
% using difference of 2 exponentials fit
%  uIPSC = exp(x./-td)-exp(x./-tr);
% INPUT:
%  D  = tr
bdebug = 0;
global target_r1090;
global td;
global dt; % Note: dt should be small enough allow sufficient accuracy in 10-90 rise estimation

% %%for debug
% % tr =.92;%% Rise  (RS basket cell data from LG on 060907)
% % td =7.66;%% Decay
tr = D(1);
% td = D(2); 
xstep = dt*1000; % ms
x = [0:xstep:80]; %% time ms
uIPSC = exp(x./-td)-exp(x./-tr);
uIPSC = uIPSC/max(uIPSC);

if bdebug
    plotdata(uIPSC',xstep/1000)
end

ind_pk = find(uIPSC == max(uIPSC));
r1090 = interp1(uIPSC(1:ind_pk),x(1:ind_pk),.9) - interp1(uIPSC(1:ind_pk),x(1:ind_pk),.1); %get 10-90 rise time
% d1090 = interp1(uIPSC(ind_pk:end),x(ind_pk:end),.1) -
% interp1(uIPSC(ind_pk:end),x(ind_pk:end),.9); %get 90-10 decay time

f = abs((target_r1090-r1090))/target_r1090; %+ (target_d1090-d1090)/target_d1090;% error function
