%% Voltage(t) 
    % INJECT Experimental currents into RC model
%% 
   global gEsyn;
global gIsyn;
global dt;
output = importStruct_abf('C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\2006_08_09_0012-1baseline.abf',1,0);
dt = output.dt;

Vdr = 90e-3;
gEsyn = output.data(:,3).*10^-12/Vdr;
gEsyn = 0.5*gEsyn; %% scale down to be more typical ration Ex/In
Vdr = 90e-3;
gIsyn = output.data(:,2).*10^-12/Vdr;
xtime = [1:length(gIsyn)]*output.dt;
%% to check RC response of cell
% gEsyn(1:(100e-3/output.dt)) = 0; %% to see response of the cell w/o currents
% gIsyn(1:(100e-3/output.dt)) = 0;

figure(1) %% CONDUCTANCE
clf
plot(xtime,gEsyn*1e9,'r')
hold on
plot(xtime,gIsyn*1e9,'b')
ylabel('Conductance (nS)')
xlabel('Time (s)')

%% assume current is conductance (driving force is the same for both Ex and
%% In)

% VINCell = -55e-3  %% USER ENTERED based on experiment
JunctP = 0e-3;
% VDriv = VINCell + JunctP;

%% INITIAL Conditions
Vi = -60e-3;

%%stepsize
tstep = dt;  %match step size to sampling rate so that there is a gsyn for each step;
tspan = [tstep:tstep:3.5]; % time in ms
y0 = [Vi];
[T,Y] = ode23(@rcNeuronExIn,tspan,y0);

figure(3)
clf
plot(T,Y*1e3,'-k')
xlabel('time (ms)')
ylabel('Membrane Potential (mV)')

hold on 
plot(xtime,gIsyn*1e9+mean(Y*1e3)-10,'-b')
plot(xtime,gEsyn*-1e9+mean(Y*1e3)-10,'-r')

%%% PLOTTING Compare EX with EX + 1nS
% figure(4)
% plot(T,Y*1e3,'-k')
% hold on 
% plot(T,Yplus*1e3,'-r')
% xlabel('time (ms)')
% ylabel('Membrane Potential (mV)')
% 
% figure(5)
% clf
% plot(T,(Yplus - Y)*1e3,'-r')
% hold on
% plot(T,(Y*1e3-mean(Y*1e3))/20 + 1.5,'-b')
