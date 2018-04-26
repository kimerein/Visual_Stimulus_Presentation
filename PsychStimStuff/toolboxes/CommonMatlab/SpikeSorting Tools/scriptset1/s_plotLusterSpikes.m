% s_plotLusterSpikes
% use with s_plotklusterSpikes or s_plotclusterspikes or
% s_plotmanualclusterspikes
%%
%%  plots spikes from cluster or kluster

clear meanSpike; clear single_unit; clear spikes_to_plot; clear unit_Nspikes; clear temp_IntraIN; clear temp_IntraEX; clear tempbINSW; clear temp_spikeData;
for k = 1:1:Nwindows % Number of figures
    figure(k+figIND)
    set(gcf,'Position',[200 0 (400+NP*200) 1120])
        orient tall
    for tt = 1:NN %% Number of subplots per figure
        j = (k-1)*NN + tt;
        if j <= MAXJ
            tsubp = (tt-1)*NP+1;
            subplot(NN,NP,tsubp:tsubp+1)
            switch (str)
                case 'manclu'
                    single_unit{j} = unit{j};
                case 'clu'
                    single_unit{j} = int32(find(T==j));
                case 'klu'
                    single_unit{j} = j;
            end
            unit_Nspikes(j) = 0;  temp = []; temp2 = temp;
            for i=1:length(single_unit{j}) % Sequentially plots spikes from each cluster in a Unit
                if i ==1
                    hold off;
                end
                temp = [temp allElectrodesXspikes(:,klusters==single(single_unit{j}(i)))];
                %% Extract spiketimes and sweep
%                 temp2 = [temp2;spikeData(klusters == single_unit{j}(i),1:3)];

            end
            
            meanSpike(:,j) = mean(temp,2);
            spikes_to_plot{j} = temp;
            unit_Nspikes(j) = unit_Nspikes(j) + size(temp,2);
            temp_spikeData{j} = temp2;
            
            plot([1:1:size(allElectrodesXspikes,1)].*Ts,spikes_to_plot{j},'k')
            axis tight
            hold on
            m= min([size(single_unit{j},1) 5]);
            temp = sprintf('Unit:%d Spikes:%d C:%s',j,unit_Nspikes(j),num2str(single_unit{j}(1:m),'%d,'));
            title(temp)


            subplot(NN,NP,tsubp+2)
            binsize = 1e-3; %(sec)
            range = .04;% range(sec) of to compute cross correlation

            title([num2str(temp) ' spikes in cluster'])
            xlabel('[sec]')
            ylabel('[mV]')
%             unit_spiketimesX = single(spikeData(klusters==j,1));
            unit_spiketimesX = temp_spikeData{j}(:,1);
            
            nSperBin =       1/ output.dt/1000 * binsize * 1000;
            xcorrBinRange = range/binsize;
            cross_corr = spiketime_xcorr(unit_spiketimesX,unit_spiketimesX,nSperBin,xcorrBinRange);
            plot(([1:size(cross_corr,2)]-range/binsize)*binsize*1000-1,cross_corr,'-g');% X axis time in Sec
            axis tight
            ylim([0 .1])

            %% ADD spike rate, EI?,
             subplot(NN,NP,tsubp+4) %% DISTANCE FROM MEAN
            Xstd = 2;
            peakmean = min(meanSpike(:,j));
            dist=sum(meanSpike(:,j)*ones(1,size(spikes_to_plot{j},2))-double(spikes_to_plot{j}),1)/size(meanSpike(:,j),1);
            s = std(dist);
            m = mean(dist);

            hold off;
            [a b]=hist(dist,100);
            plot(b,a,'-k');
            hold all;
                line([Xstd*s Xstd*s],[0 max(a)]);
                line([-Xstd*s -Xstd*s],[0 max(a)]);
%                 line([-peakmean -peakmean],[0 max(a)],'Color','r');
                Dirty(j) = 2*s/abs(peakmean);
            title(['distance from mean   Dirty:' num2str(Dirty(j),'%1.2f')]);
            hold off;
            axis tight
            %% PLOT average Intracellular signal
            subplot(NN,NP,tsubp+3)
            kk =1;
            if nIntracellsignals >0
                % correct index
                    %%plots only first 1000 spikes to save memory)
                     temp1 = zeros(size(temp_spikeData{j},1),int32(Ts*60e-3)+1,nIntracellsignals); %% for EX
                  temp1 = zeros(size(temp_spikeData{j},1),int32(Ts*60e-3)+1,1); %% for EX
                    i0 = 1; i1 = 1;
                    for jjjj=1:(size(temp_spikeData{j},1))
                        temp2 = (temp_spikeData{j}(jjjj,1));
                        %% find if EX or IN
                        sweepnum = floor(double(temp2)/(N_samples/N_chn)) + 1;
                        if ~isempty(lastINSW)
                            if  sweepnum > lastINSW %% if 0mV is first in file
                                tempbINSW{j}(jjjj)=0;
                            else
                                tempbINSW{j}(jjjj)=1;
                            end
                        else
                            if  sweepnum > lastENSW %% if 0mV is last in file
                                tempbINSW{j}(jjjj)=1;
                            else
                                tempbINSW{j}(jjjj)=0;
                            end
                        end
                        % catch cases where intradata doesn't
                        % exist(beggining or end of file)
                        if (temp2-int32(Ts*10e-3)) <1  % BEGINNING
                            spikeIntraDataMissing(kk) = temp2;
                            kk = kk+1;
                            temp1(jjjj,:) =  [IntraCellData(1)*ones(1,abs(temp2-int32(Ts*10e-3))+1,'int16') IntraCellData((1):(temp2 + int32(Ts*50e-3)))];
                        elseif (temp2 + int32(Ts*50e-3)) > size(IntraCellData,1)*size(IntraCellData,2) % END
                            spikeIntraDataMissing(kk) = temp2;
                            kk = kk+1;
                            temp1(jjjj,:,1) =  [IntraCellData((temp2-int32(Ts*10e-3)):end) [IntraCellData(end)*ones(1,(temp2+int32(Ts*50e-3))-size(IntraCellData,1)*size(IntraCellData,2),'int16')]];
                        else
                            temp1(jjjj,:,1) = IntraCellData((temp2-int32(Ts*10e-3)):(temp2 + int32(Ts*50e-3)));
                        end
 
                    end
                    IntraCHN = 1;
                        temp_Intra{j,IntraCHN} = single(temp1)*single(output.gain(IntraCHN))+ single(output.offset(IntraCHN));
  
                    %                     plot([1:size(temp1(jjjj,:),2)]*1e3/Ts,temp_Intra{j},'
                    %                     -k'); %% plot every intrasignal
                    %                        plot([1:size(temp1(jjjj,:),2)]*1e3/Ts,(mean(temp1,1)+std(single(temp1),1))*output.gain(jjj),'--k');
                    %                                        hold on;
                    hold off
                    [r c]= find(tempbINSW{j}==1);
                    temp = mean(temp_Intra{j}(c,:),1) - mean(mean(temp_Intra{j}(c,:),1)) + 10;
                    plot([1:size(temp1(jjjj,:),2)]*1e3/Ts,temp,'-r');
                    hold on
                    [r c]= find(tempbINSW{j}==0);
                    temp = mean(temp_Intra{j}(c,:),1) - mean(mean(temp_Intra{j}(c,:),1));
                    plot([1:size(temp1(jjjj,:),2)]*1e3/Ts,temp,'-b');
                    axis tight
                    %                     ylim([-50+min(temp) 50+max(temp)]);
                    xlabel('ms')
                    ylabel('pA')
                    temp = tempbINSW{j};
                    title(['INr: ' num2str(sum(temp)) ' EXr: ' num2str(size(temp,2)-sum(temp))])
              
                if b_intraspikes
                    subplot(NN,NP,tsubp+4)
                    binsize = 5e-3;
                    range=0.1;
                    %         nbins =        (N_samples/N_chn * sw)/
                    %         Ts/(binsize);
                    nSperBin =       1/ output.dt/1000 * binsize * 1000;
                    xcorrBinRange = range/binsize;
                    unit_spiketimesX = single(spikeData(klusters==j,1));
                    cross_corr = spiketime_xcorr(first_InspikeData(:,1),unit_spiketimesX,nSperBin,xcorrBinRange);
                    plot(([1:size(cross_corr,2)]-range/binsize)*binsize*1000-1,cross_corr,'-b');% X axis time in Sec
                    axis tight
                end
            end

        end
    end
     if breport
         dirtemp = 'REPORT';
         figdesc = [figString '_EE' num2str(triggerEE) '_W' num2str(k)];
         savefigure(writedirheader,dirtemp,figdesc,'png')
     end


end

if breport
    %%save data to file
    dirtemp = '';
    spath = [writedirheader dirtemp '\'];
    if (~isdir(spath))
        mkdir(spath);
    end
    temp = sprintf('%s%sEE%d',spath,figdesc,triggerEE);
    %     print('-dtiffn',temp);
    save(temp,'sourcefile','EEgroup','klusters','Center','num_c','temp_spikeData','temp_Intra','tempbINSW','spikes_to_plot','unit_Nspikes','std_threshold','minTime','Ts')  %% NOTE font size too big on axis (overlaps)
end
