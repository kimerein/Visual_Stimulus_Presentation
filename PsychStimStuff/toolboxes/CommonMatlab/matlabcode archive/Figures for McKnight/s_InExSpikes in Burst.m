%%% Plotting Intra and EE spikes from a Burst for McKnight

nIntracellsignals =1;
indexoffset =1;
std_threshold = 6;
maxpeakamp = 10000*std_threshold;  %in std
maxpeakwidth =  2e-3;
WOI = [750e-6 750e-6; 1500e-6 1500e-6];%Window of interest around Beg_spike -1000us + 2000us

%% LOAD data from 
%% Cell 7/14/2005 0017  PYR with extracellular burst
output = importStruct_abf('C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\Rotation 2005 data\2005_07_14_0017.abf',1,1)
%% Just one sweep extracted

intradata = output.data(:,indexoffset+nIntracellsignals);
dintradata = - diff(intradata);
eedata = output.data(:,indexoffset+nIntracellsignals +1);
intradata = intradata(2:end);
eedata = eedata(2:end);
% Extract Intra Spikes
Intrastd(nIntracellsignals) =  mean(std(single(output.data(:,indexoffset+nIntracellsignals:output.Nchan:end))));
Intramean(nIntracellsignals) =  mean(mean(single(output.data(:,indexoffset+nIntracellsignals:output.Nchan:end))));
temp = int32(FindSpikes(intradata, Intramean(nIntracellsignals)+std_threshold*Intrastd(nIntracellsignals) ,1/output.dt,1));

Intra = zeros(size(temp,1),251); EE = Intra;dVdt = EE;
try
for i =1:size(temp,1)  % extract spikes
  temp1 = intradata((temp(i)):temp(i)+150);
  temp2 = find(temp1 == max(temp1));
  Intra(i,:) = intradata(temp(i)+temp2-100:temp(i)+temp2+150);
   EE(i,:) = eedata(temp(i)+temp2-100:temp(i)+temp2+150);
   dVdt(i,:) = dintradata(temp(i)+temp2-100:temp(i)+temp2+150);
end
catch
    i
end

%% normalized spikes, 1-X in a burst
Intra = Intra'./max(max(Intra));
EE = -EE'./min(min(EE));
dVdt = -dVdt'./min(min(dVdt));

figure(1);
x = [1:size(dVdt,1)].*output.dt*1000;
i =1; %% plot first spike
plot(x,Intra(:,i),'-k',x,EE(:,i),'-b',x,dVdt(:,i),'-r');
title('First Spike');
axis tight
figure(2);
i = size(Intra,2); %% last spike
plot(x,Intra(:,i),'-k',x,EE(:,i),'-b',x,dVdt(:,i),'-r')
title('Last Spike');
ylim([-1 1])
xlim([0 5])

%% Plot whole trace
figure(3);
x = [1:size(intradata,1)].*output.dt*1000;
plot(x,intradata);
axis tight
figure(4)
plot(x,eedata);
axis tight

