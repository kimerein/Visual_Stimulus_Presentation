s_loadSynEvents


% Index of event to be analyzed
A=7;        %% USER ENTRY
B=8;        %% USER ENTRY
breport = 1;

defineDir
expttype = 'KinateOsc';
writedirpath = [DATAANAL_DIR expttype '\'];
writedirheader = [DATAANAL_DIR expttype '\'];

%% RUN BY HAND
s_loadSynEvents
% debug
% a = [1:10];
% [event(A).time(ind(a,1)),event(A).peak(ind(a,1)),event(B).time(ind(a,2)),event(B).peak(ind(a,2))]
% plot(event(A).peak(ind(:,1)),event(B).peak(ind(:,2)),'.b')

%% check whether peak amplitude predicts next peak amplitude
s_SynEvent_ExInScale


figure(2)
%% not in real units
[a x] =hist(diff(event(A).peak),100)
 stairs(a,'r','LineWidth',2);
 hold on
 [a x] =hist(diff(event(B).peak),100)
 stairs(a,'b','LineWidth',2);
 title('Histogram delta event Peak')

s_SynEventsStats