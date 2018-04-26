function param = xSpikeParamsCon(data,dt,nElectrodes)
% function param = xSpikeParamsCon(data,dt,nElectrodes)
% same as xSpikeParams but with concatanated spikes
% dt = time between samples 
% nElectrodes in concatanated data.
% Spike Parameters
% amplitude of P1,T1,P2
% width of P1,T1,P2
%% assume DC component is removed.
%% assume spikes take form Peak(P1),Trough(T1),Peak(P2).
%% parm = [P1 T1 P2 P1W T1W P2W T1F T1R P2R P2F];
nPoints = size(data,1)/nElectrodes;
for i=1:nElectrodes
    param{i} = xSpikeParams(data((i-1)*nPoints+1:i*nPoints,:)',dt,0);
end