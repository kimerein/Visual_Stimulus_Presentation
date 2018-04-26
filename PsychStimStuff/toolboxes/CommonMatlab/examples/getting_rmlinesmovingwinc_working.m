    %     remove 60Hz
    %     df(:,ch_LFP,:) = prepdata(squeeze(d(:,ch_LFP,:))',0,dt,[NaN NaN 1])*1000; % for 0325025
    %     params.tapers = [3 5];params.Fs  = 1/dt;params.fpass = [0 500];
    %     movingwin = [10 50]; %sec
    %     tau = 10; p = 0.05;f0=[];
    % [datac,datafit,Amps,freqs]=rmlinesmovingwinc(squeeze(d(:,ch_LFP,:))',f0,movingwin,tau,params,p,'y');
