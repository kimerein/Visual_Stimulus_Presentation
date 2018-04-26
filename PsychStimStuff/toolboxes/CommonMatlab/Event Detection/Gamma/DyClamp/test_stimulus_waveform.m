%% test conductance injection

output = importStruct_abf('C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\2006_10_22_0013.abf',-1,0);
I = output.data(2216:2216+200,3);
V_sig = abs(output.data(2216:2216+200,4));
figure(1)
clf
plot([1:size(I)]*output.dt*1000,I/max(I),'-b')
hold on
plot([1:size(V_sig)]*output.dt*1000,V_sig/max(V_sig),'-r')
xlabel('ms')

I = output.data(4438:4438+200,3);
V_sig = -1*output.data(4438:4438+200,4);
plot([1:size(I)]*output.dt*1000,I/max(I),'-g')

plot([1:size(V_sig)]*output.dt*1000,V_sig/max(V_sig),'-k')
