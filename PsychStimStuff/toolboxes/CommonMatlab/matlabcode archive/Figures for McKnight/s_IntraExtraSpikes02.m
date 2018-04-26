
%plot intra, extra and dvdt for McKnight grant
% output = importStruct_abf('C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\subset2006_03_20_0009.abf',-1,1);
% PYR
% output = importStruct_abf('C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\2005_09_28_0056.abf',-1,1);
%FS
tic
output = importStruct_abf('C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\2005_12_20_0005.abf',-1,1);
toc
N_chn = output.Nchan;
Ts = 1/output.dt;

% define window where intracellular spikes occur
data = output.data(31430:31710,:); % Determined expt by expt
temp1 = data(:,indexoffset+nIntracellsignals:N_chn:end); 
% exclude sweeps
temp = [1:15 17:42 44:51 53:61 63:98];
intradata = temp1(:,temp);
temp1 = data(:,indexoffset+nIntracellsignals+3:N_chn:end); 
eedata = temp1(:,temp);

% Extract Intra Spikes
nIntracellsignals =1;
indexoffset =5;
std_threshold = 6;
maxpeakamp = 10000*std_threshold;  %in std
maxpeakwidth =  2e-3;
WOI = [750e-6 750e-6; 1500e-6 1500e-6];%Window of interest around Beg_spike -1000us + 2000us
Intrastd(nIntracellsignals) =  mean(std(single(output.data(:,indexoffset+nIntracellsignals:N_chn:end))));
Intramean(nIntracellsignals) =  mean(mean(single(output.data(:,indexoffset+nIntracellsignals:N_chn:end))));
temp = int32(FindSpikes(intradata, Intramean(nIntracellsignals)+std_threshold*Intrastd(nIntracellsignals) ,Ts,1));

Intra = zeros(size(temp,1),61); EE = zeros(size(temp,1),61);
try
for i =1:size(temp,1)  % extract spikes
  temp1 = intradata((temp(i)-20):temp(i)+40);
  temp2 = find(temp1 == max(temp1));
  Intra(i,:) = intradata(temp(i)+temp2-40:temp(i)+temp2+20);
   EE(i,:) = eedata(temp(i)+temp2-40:temp(i)+temp2+20);
end
catch
    i
end
EE = EE'; Intra = Intra';
% EE  = E(:,[1:10,12:end]);
% EE = EE(480:780,:);
mEE = mean(EE,2);
% subtract baseline
mEE = mEE - mean(mEE(1:10),1);

% Intra = Intra(:,[1:10,12:end]);
% Intra = Intra(480:780,:);
mI = mean(Intra,2);
x = [2:size(mI,1)]';
dVdt = (-diff(mI));
dVdt = dVdt - mean(dVdt(1:10),1);
dVdt = dVdt/min(dVdt);
% plot(x,mI(1:end-101),'-b',x,dVdt,'-r')
figure;
plot(x,mI(2:end)/max(mI),'-k')
hold on;
plot(x,-mEE(2:end)/min(mEE),'-b',x,-dVdt(1:end),'-r')

%% for checking when the current step occured
% IC = output.data(3500:4500,7);
% IC = IC(480:780,:);

