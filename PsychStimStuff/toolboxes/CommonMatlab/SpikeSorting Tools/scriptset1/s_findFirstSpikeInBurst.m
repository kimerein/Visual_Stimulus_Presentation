spikeselection = 'Xintra';
Xstd = 1; %% number of std from mean spike distance to reject
XIstd =1.6;

NN = 4;
NP =6;

IMAX = size(unit,2);
% close(figIND);
Nwindows  = ceil(size(selectUnits,2)/NN);

for k = 1:Nwindows
    figure(figIND + k)
    set(gcf,'Position',[200 0 (400+NP*200) 1120])
        orient tall
    for subi = 1:NN
        ii = subi + NN*(k-1);
        if ii < IMAX+1
            i = selectUnits(ii);
            %         if bXstragglerspikes
            unit(i).dist=sum(mean(unit(i).spikes,2)*ones(1,size(unit(i).spikes,2))-double(unit(i).spikes),1);
            s = std(unit(i).dist);
            m = mean(unit(i).dist);
            spike_ind = find(abs(unit(i).dist)<abs(m-Xstd*s)) ;%find spikes within Xstd of mean
 
            unit(i).Intradist=sum(mean(unit(i).intraData,2)*ones(1,size(unit(i).intraData,2))-double(unit(i).intraData),1);
            s = std(unit(i).Intradist);
            m = mean(unit(i).Intradist);
            spikeIntra_ind = find(abs(unit(i).dist)<abs(m-XIstd*s)) ;%find spikes within Xstd of mean

            % Include only spikes within Xstd of mean
            tempIND = [1:size(unit(i).spikes,2)];
            stemp = ['all Spikes'];
            switch (spikeselection)
                case 'first'
                tempIND = find(unit(i).spikeOrder==1);
                stemp = ['1stSpk'];
                case 'Xintra'
                 tempIND = intersect(intersect(find(unit(i).spikeOrder==1),spike_ind),spikeIntra_ind);
                stemp = ['< ' num2str(Xstd) ' std&1stSpk' ' &' num2str(XIstd) '< Intra'];   
                case 'Xspk'
                tempIND = intersect(find(unit(i).spikeOrder==1),spike_ind);
                stemp = ['< ' num2str(Xstd) ' std&1stSpk'];
            end
            subplot(NN,NP,1+(subi-1)*NP) %% EE WAVEFORM

            xtime = 1000*dt.*[1:size(unit(i).spikes,1)];
            hold off;
            plot(xtime, unit(i).spikes(:,tempIND),'-k')
            hold on
            mspike=mean(unit(i).spikes(:,tempIND),2)' ;
            plot(xtime,mspike,'-r')

            %% get FWHM Peak2 from electrode with largest mean spike
            ind_min = find(mspike==min(mspike));
            electrode_min = floor(ind_min/(size(mspike,2)/size(EEgroup,2))) +1;
            unit(i).SpikeParam = xSpikeParams2(mspike((electrode_min-1)*(size(mspike,2)/size(EEgroup,2))+1:(electrode_min)*(size(mspike,2)/size(EEgroup,2))-1),1/Ts*1000);
            spikeParam(i) = unit(i).SpikeParam(12) - unit(i).SpikeParam(11);
            title(['AUnit: ' num2str(i) ' P2W:' num2str(spikeParam(i)) ' N:' num2str(sum(unit(i).spikeOrder)) '(' num2str(size(spike_ind,2)) ')'])
            axis tight
            if subi==NN
                xlabel('ms')
            end

            figure(figIND + k) %% incase it got changed somewhere

%             subplot(NN,NP,2+(subi-1)*NP) %% XCORR
% 
%             binsize = 1e-3; %(sec)
%             range = .04;% range(sec) of to compute cross correlation
%             hold off;
%             unit_spiketimesX = single(unit(i).spikeData(tempIND));
%             nSperBin =       Ts * binsize ;
%             xcorrBinRange = range/binsize;
%             cross_corr = spiketime_xcorr(unit_spiketimesX,unit_spiketimesX,nSperBin,xcorrBinRange);
%             plot(([1:size(cross_corr,2)]-range/binsize)*binsize*1000-1,cross_corr,'-g');% X axis time in Sec
%             axis tight
%             ylim([0 .1])
%             xlabel('[msec]')
%             ylabel('[mV]')
%             title('Xcorr')

            subplot(NN,NP,2+(subi-1)*NP) %% JITTER
            hold off;
            bIntra =1;
            bIPSC = 0; %% Excitatory
            [r c1] = find(unit(i).bIN==1);
            temp = intersect(tempIND,c1);
            [PSP_time Failure] = PSPParam(unit(i).spikeData(temp,1),unit(i).intraData(:,temp),1/Ts,bIPSC,bIntra);
            [b a] = hist(PSP_time(:,2),80);
            bar(a,b,'r');
            temp1 = max(b(2:end));
            hold on
            bIPSC = 1; %% Inhibitory
            [r c0] = find(unit(i).bIN==0);
            temp = intersect(tempIND,c0);
            [PSP_time Failure] = PSPParam(unit(i).spikeData(temp,1),unit(i).intraData(:,temp),1/Ts,bIPSC,bIntra);
            [b a] = hist(PSP_time(:,2),80);
            bar(a,-b,'b');
            axis tight;
            xlim([0 10])
            ylim([-max(b(2:end)) temp1]);
            title('PsP Jitter')
            if subi==NN
                xlabel('ms')
            end
            subplot(NN,NP,3+(subi-1)*NP)  %% INTRACELL WAVEform

            bIPSC = 0;
            xtime = dt*1000.*[1:size(unit(i).intraData,1)]-10;
            temp = mean(unit(i).intraData(:,intersect(tempIND,c1)),2);
            stemp1 = PSPParam(1,temp,1/Ts,bIPSC,bIntra);
            hold off
            plot(xtime, temp,'-r')
            axis tight
            bIPSC = 1;
            xtime = dt*1000.*[1:size(unit(i).intraData,1)]-10;
            temp = mean(unit(i).intraData(:,intersect(tempIND,c0)),2);
            stemp0 = PSPParam(1,temp,1/Ts,bIPSC,bIntra);
            hold on
            plot(xtime, temp,'-b')
            %             plot(xtime, mean(unit(i).intraData(:,setdiff(Failure,itemp)),2),'b');
            axis tight
            %             ylim([-50+min(temp) 50+max(temp)]);
            if subi==NN
                ylabel('pA')
                xlabel('ms')
            end
            %%% Calculate Latency from Mean
            title([stemp ' 2PkInh:' num2str(stemp0(:,2),'%1.1f') ' 2PkEx:' num2str(stemp1(:,2),'%1.1f')]);

            subplot(NN,NP,4+(subi-1)*NP) %% SPIKE RATE
            hold off;
            [r c] = find(unit(i).spikeData(:,1)~=0);
            temp = single(unit(i).spikeData(intersect(tempIND,r),1));
            [b a] = hist(temp,max(temp)/(Ts)) ;%% spike rate over time
            plot(a./60./Ts,b);
            if subi==NN
                ylabel('spikes/sec')
                xlabel('Mins (1sec bins)')
            end
            title('Spike Rate');
            hold off;
            axis tight
            
            subplot(NN,NP,5+(subi-1)*NP) %%ISI
            hold off;    
            temp = diff(single(unit(i).spikeData(intersect(tempIND,r),1)));
            [a b] = hist(temp,max(temp)/(Ts/1000)); %% ISI
            plot(b,a,'-k');
            if subi==NN
                ylabel('ms')
                xlabel('1ms bins')
            end
            title('ISI');
           axis tight
 
            
            subplot(NN,NP,6+(subi-1)*NP) %% DISTANCE FROM MEAN SPIKE
            hold off;
            [a b]=hist(unit(i).dist,100);
            plot(b,a,'-k');
            hold all;
            if bXstragglerspikes
                line([Xstd*s Xstd*s],[0 max(a)]);
                line([-Xstd*s -Xstd*s],[0 max(a)]);
            end
            title('distance from mean');
            hold off;
            axis tight

        end
    end
    if breport
        dirtemp = 'REPORT';
        figdesc = [figString '_EE' num2str(triggerEE) '_W' num2str(k)];
        savefigure(writedirheader,dirtemp,figdesc,'png')
    end

end

figure(570)
hist(spikeParam)
title('Histogram: Time From Trough to Peak')
xlabel('ms')

if breport
    figString = 'ASpikes_TtoP_Hist';

    dirtemp = 'REPORT';
    figdesc = [figString '_EE' num2str(triggerEE)];
    savefigure(writedirheader,dirtemp,figdesc,'png')
end

figure(571)
hist(Dirty,20)
title('KlusterDirtyHist')


% if bsave
%         save(temp,'klusters','Center','num_c')  %% NOTE font size too big on axis (overlaps)
