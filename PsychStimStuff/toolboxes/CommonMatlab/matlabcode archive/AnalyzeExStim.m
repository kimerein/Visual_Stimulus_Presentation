%***********************************************************************
% AnalyzeExStim.m
%
% Flavio Frohlich
% Last change: 7/10/2005
% Reads .atf file, basic analysis of extracellular traces with activity
% evoked by extracellular stimulation.
%
%***********************************************************************




% Initialize workspace

close all;
clc;

my_spikes = [];
candidates_proc = [];
my_ind = [];

Ts =  1/49e-6; %1/20e-6; %% Sampling Rate

filenumber = 26;
%my_data = dlmread(['SCZ\stim\e057080' num2str(filenumber) '.atf'],'\t',12,0);
% my_data = dlmread(['057080' num2str(filenumber) '.atf'],'\t',12,0);

filenumber_arte = 28;
% my_data_arte = dlmread(['057080' num2str(filenumber_arte) '.atf'],'\t',12,0);
[duration N] = size(my_data);



% -------------------------------------------------------------------------
% Plot traces

my_beg =  0.095; %0.5015;
my_end =  0.110; %0.530;

chn = [2:N]; % Note: Chn 1 is time

figure(1)
plot(my_data(my_beg*Ts+1:my_end*Ts,1),my_data(my_beg*Ts+1:my_end*Ts,chn))
ylabel('mV')
ylim([-0.02 0.02])





% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% Compute std
stdtrial = std(my_data(1:0.1*Ts,2:end));
stdall = mean(stdtrial)




ttx = 1

if ttx ==1

% -------------------------------------------------------------------------
% Subtract averaged stimuls artefact (whith bath applied TTX)

my_beg = 0.095;
my_end = 0.110;

chn = [2:2]; % Note: Chn 1 is time

my_data_proc = my_data;
artefact = mean(my_data(:,[2:N]),2);

figure(2)
plot(my_data(my_beg*Ts+1:my_end*Ts,1),artefact(my_beg*Ts+1:my_end*Ts))


my_data_proc(:,2:N) = my_data_proc(:,2:N)-artefact*ones(1,N-1);

chn = [2:N]
figure(66)
plot(my_data(my_beg*Ts+1:my_end*Ts,1),my_data_proc(my_beg*Ts+1:my_end*Ts,chn))
ylabel('mV')
ylim([-0.02 0.02])


figure(3)
title('Subtracted Stim Artefact')
plot(my_data(my_beg*Ts+1:my_end*Ts,1),my_data_proc(my_beg*Ts+1:my_end*Ts,chn))
ylabel('mV')
hold on
% -------------------------------------------------------------------------

end

slpa = 0;
if slpa == 1

% -------------------------------------------------------------------------
% Use SALPA to remove stim artefact

% Initialize variables

a0 = []
Tk = [];
ttemp = 0;
W = [];

NSamples = size(my_data,1);
NSalpa = 10;%50;
Vall = zeros(NSamples-2*NSalpa,N);
start = 1;%2080; %round(0.500*Ts);% 2080; % after depegging
NSamples = NSamples -start + 1;

temp = [-NSalpa:1:NSalpa];

for i=0:1:6
    Tk = [Tk sum(temp.^i)];
end

S = [Tk(1) Tk(2) Tk(3) Tk(4); Tk(2) Tk(3) Tk(4) Tk(5);Tk(3) Tk(4) Tk(5) Tk(6);Tk(4) Tk(5) Tk(6) Tk(7)];
S=inv(S);

for sweep = 2:1:2
    V = my_data(start:end,sweep)';


    for j=NSalpa+1:1:NSamples-2*NSalpa;

        for ind=0:1:3
            W(ind+1)= temp.^ind*V(j-NSalpa:j+NSalpa)';
        end
        a0(j-NSalpa) = S(1,1:end)*W';

        if j==NSalpa+1
            a0_beg = a0(1);
            a1_beg = S(2,1:end)*W';
            a2_beg = S(3,1:end)*W';
            a3_beg = S(4,1:end)*W';
            Wold = W
        end
    end


    Vnew_beg = [];
    for i=1:NSalpa
        Vnew_beg(i) = V(i) - (a0_beg +a1_beg*(i-NSalpa-1)+ a2_beg*(i-NSalpa-1)^2 + a3_beg*(i-NSalpa-1)^3);
    end

    Vnew = V(NSalpa+1:end-2*NSalpa) - a0;

    Vall(:,sweep) = [my_data(1:start-1,sweep)' Vnew_beg Vnew]';

    figure(101)
    plot(Vall)
    figure(3)
    hold on
    plot(my_data(my_beg*Ts+1:my_end*Ts,1),Vall(my_beg*Ts+1:my_end*Ts,sweep),'r')
    % line([my_data(round(my_beg*Ts),1) my_data(round(my_end*Ts),1)],[stdall stdall].*(-3),'Color','k')
    % line([my_data(round(my_beg*Ts),1) my_data(round(my_end*Ts),1)],[stdall stdall].*(-5),'Color','k')
       plot(my_data(my_beg*Ts+1:my_end*Ts,1),my_data(my_beg*Ts+1:my_end*Ts,sweep),'k')
    % ylim([-5e-3 4e-3])
    if ttx == 1
    legend('TTX Subtraction','SALPA')
    end
end
end
% -------------------------------------------------------------------------


spike_sort = 1


if spike_sort == 1

% -------------------------------------------------------------------------
% Find and extract spikes

wof_beg =  0.105;% 0.502;
wof_end =  0.1075; %0.530;


spikecount = zeros(N-1,1);
for i=2:1:N
    candidates_proc = [];
    neg_tr = -2e-3;
    clear candidates;
    % candidates = find((my_data_proc(:,i) > pos_tr));
    
    % my_data_proc = Vall;
    candidates = find((my_data_proc(:,i) < neg_tr))
    
    candidates = candidates(candidates>wof_beg*Ts);
    candidates = candidates(candidates<wof_end*Ts);
    if (isempty(candidates) ~=1)
        candidates_proc = candidates(1);
        for j=2:length(candidates)
             if candidates(j)-candidates(j-1) > (Ts*0.002)
            
                candidates_proc = [candidates_proc; candidates(j)];
                
             end
        end




        for k=1:length(candidates_proc)
            if candidates_proc(k) ~= -1
                [dummy ind] = min(my_data_proc(candidates_proc(k)-(Ts*0.0005):(candidates_proc(k)+(Ts*0.003)),i));
              
                % align on treshold crossing
                % my_spikes = [my_spikes my_data(candidates_proc(k)-(Ts*0.0005):(candidates_proc(k)+(Ts*0.002)),i)];
                % align on negative peak
                
                offset = 20;
                  my_spikes = [my_spikes my_data_proc(candidates_proc(k)-(Ts*0.0005)+ind-offset:(candidates_proc(k)+(Ts*0.002))+ind,i)];
                  my_ind = [my_ind (candidates_proc(k)+ind)/Ts - 0.10055];
                  spikecount(i-1) =  (candidates_proc(k)+ind)/Ts - 0.10055;
                
            end
        end
    end
end
figure(4)
plot(my_spikes)
figure(5)
hist(my_ind)
    

[PC, SCORE, LATENT, TSQUARE] = princomp(my_spikes');
figure(50)
plotmatrix(SCORE(:,1:5))
cluster1 = [];
cluster2 = [];

clusters = kmeans(my_spikes',2)

cluster1_ind = find(clusters == 1);
cluster2_ind = find(clusters == 2);

cluster1 = my_spikes(:, cluster1_ind);
cluster2 = my_spikes(:, cluster2_ind);

% 
figure(51)

subplot(2,1,1)
plot([1:1:length(cluster1(:,1))]./Ts,cluster1)
ylabel('[mV]')

subplot(2,1,2)
plot([1:1:length(cluster2(:,1))]./Ts,cluster2)
ylabel('[mV]')
xlabel('Time [s]')

figure(52)
plotmatrix(SCORE(cluster1_ind,1:5))
figure(53)
plotmatrix(SCORE(cluster2_ind,1:5))

% % 
% % % -------------------------------------------------------------------------
% % 
% % 
% % 
end
