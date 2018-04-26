%% script to read manually entered synaptic event data
%% into structure 
%% USER DEFINE
% spsh according to:
%
%% spsh has columns:
%Trace	Event Start Time (ms)	Rise Tau (ms)	Decay Tau (ms)	Peak Amp (pA)	Interevent Interval (ms)

spsh = single(spsh);                            %% USER DEFINE
filename = '2006_08_09_0012ss_InEventIn1';        %% USER ENTRY
V = -3;                                         %% USER ENTRY
expnum = 's1';                                  %% USER ENTRY
savefilename = ['2006_08_09_' expnum];          %% USER ENTRY

if exist('event')
i = size(event,2)+1;
else
    i = 1;
end

event(i).sfilename = filename;
event(i).V = V;
event(i).expnum = expnum;
event(i).Trace = spsh(:,1);
event(i).time = spsh(:,2);
event(i).RiseTau = spsh(:,3);
event(i).DecayTau = spsh(:,4);
event(i).peak = spsh(:,5);
event(i).isi = spsh(:,6);

%
save([writedirpath 'Pairs\' savefilename],'event')
% 
% event(i).sfilename = filename;
% event(i).V = V;
% event(i).expnum = expnum;
% event(i).Trace = event(i).Trace(1:1004);
% event(i).time = event(i).time(1:1004);
% event(i).RiseTau = event(i).RiseTau(1:1004);
% event(i).DecayTau = event(i).DecayTau(1:1004);
% event(i).peak = event(i).peak(1:1004);
% event(i).isi = event(i).isi (1:1004);
% 

