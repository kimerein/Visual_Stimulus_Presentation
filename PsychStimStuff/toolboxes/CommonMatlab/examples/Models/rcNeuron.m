function dVdt= rcNeuron(t,V)
%function dVdt = rcNeuron(V)
%% RC-Model of neuron
% 
%% Parameters: USER DEFINED BASED on experiment
R = 188.15e6;      % (Ohm) Resistance of Neuron
C = 103.89e-12;    % (F) capacitance of Neuron
%%OR
% tau = R*C;
global gsyn;
global dt;

JunctP = -5e-3;
Vm = -63.5e-3 + JunctP;    % Reversal potential of Leak Current (Resting Potential)


% %% SYNAPTIC
Stau = 10e-3 ;  % (ms);
VRse = 0 ;     % Reversal potential of Excitatory currents

Gse = 1;   % (S) max synaptic conductance


dVdt = -1/C* (1/R*(V-Vm) + gsyn(round(t/(dt)))*(V-VRse));
%
% dVdt = -1/tau*(V-Vm) -1/C*( Gse*gsyn*(V-VRse));


% 
