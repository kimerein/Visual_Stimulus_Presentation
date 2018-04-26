% s_MSPairFigure
%  selectUnits = [good_type];
spikeselection = 'Xintra';
% spikeselection = 'Xspk';
Xstd = 1; %% number of std from mean spike distance to reject
XIstd =2;
sGrid  = 5; %ms in grid
            bplotraw = 0;


Xstd =1; %% number of std from mean spike distance to reject
figIND = 800;
IMAX = size(selectUnits,2);
NU = 2;  %% number of units per window
NSU = [2 1]; %%subplots per units 
NS = NU*NSU(1);

eW = 0.08; eH = 0.035; %%edges
bH = 0.0005;% between rows
sW = 1-1.4*eW; sL =eW;
sH = 1/NS-1*eH;

NN = min(NU,IMAX); %% cells per window

figString = ['MSFig_' sunitfilename];
% close(figIND);
H = 1075; W = H*.77; % window size
Nwindows  = ceil(size(selectUnits,2)/NN);
for k = 1:Nwindows
    figure(figIND + k) %% Pair EE and Intra WINDOW
    clf
    set(gcf,'Position',[775 50 W H],'PaperUnits','normalized')
    orient tall
    for subi = 1:NN
        ii = subi + NN*(k-1);
        if ii < IMAX+1
            i = selectUnits(ii);
            spikeIntra_ind = [];
             unit(i).dist=sum(mean(unit(i).spikes,2)*ones(1,size(unit(i).spikes,2))-double(unit(i).spikes),1);
            s = std(unit(i).dist);
            m = mean(unit(i).dist);
            spike_ind = find(abs(unit(i).dist)<abs(m-Xstd*s)) ;%find spikes within Xstd of mean
 
            %% find distribution of EPSC/IPSCs around the mean (so that
            %% later the stragglers can be excluded from plots)
            [r c0] = find(unit(i).bIN==0); %% first at 1 holding potential 
            unit(i).Intradist(:,c0)=sum(mean(unit(i).intraData(:,c0),2)*ones(1,size(unit(i).intraData(:,c0),2))-double(unit(i).intraData(:,c0)),1);
            s = std(unit(i).Intradist(:,c0));
            m = mean(unit(i).Intradist(:,c0));
            spikeIntra_ind = min(c0) -1 + find(abs(unit(i).Intradist(:,c0))<abs(m-XIstd*s)) ;%find spikes within Xstd of mean
            [r c1] = find(unit(i).bIN==1);   %% now the other holding potential
            unit(i).Intradist(:,c1)=sum(mean(unit(i).intraData(:,c1),2)*ones(1,size(unit(i).intraData(:,c1),2))-double(unit(i).intraData(:,c1)),1);
            s = std(unit(i).Intradist(:,c1));
            m = mean(unit(i).Intradist(:,c1));
            spikeIntra_ind = [(min(c1) -1 +find(abs(unit(i).Intradist(:,c1))<abs(m-XIstd*s))) spikeIntra_ind] ;%find spikes within Xstd of mean
            
            % Include only spikes within Xstd of mean
            tempIND = [1:size(unit(i).spikes,2)];
            stemp = ['all Spikes'];
            switch (spikeselection)
                case 'first'
                tempIND = find(unit(i).spikeOrder==1);
                stemp = ['1stSpk'];
                case 'Xintra'
                 tempIND = intersect(intersect(find(unit(i).spikeOrder==1),spike_ind),spikeIntra_ind);
                stemp = ['< ' num2str(Xstd) ' std&1stSpk' ' & < ' num2str(XIstd) 'Intra'];   
                case 'Xspk'
                tempIND = intersect(find(unit(i).spikeOrder==1),spike_ind);
                stemp = ['< ' num2str(Xstd) ' std&1stSpk'];
            end
            
            %% calculate position of subplot
            nSU =2+(subi-1)*NSU(1); % subplot in unit
%             sB = 1-((nSU-1)/(NS+NS*bH+eH*1.2)+(nSU-1)*bH +eH);
           sB = 1-((nSU)/(NS+NS*bH+eH*5)+(nSU-1)*bH );


            subplot('Position',[sL sB sW sH]) %% INTRA
            clear mIntra;
%             xrange = [int32(5e-3*Ts):int32(20e-3*Ts)];% s
            xrange = [int32(5e-3*Ts):int32(60e-3*Ts)];% s
            xtime1 = dt*1000.*single(xrange)- 10;
            bIntra =1;
            %% IN
            bIPSC = 1;  
            shiftmIntra0 = 10;
             stemp0 = [0 0 0]; mIntra0 = shiftmIntra0;
            temp00 = intersect(tempIND,c0);
            if ~isempty(temp00)
                mIntra0 = mean(unit(i).intraData(:,temp00),2);
%                 sIntra0 =  std(unit(i).intraData(:,temp00),0,2);
 
                temp0 = mean(mIntra0(1:int32(10e-3*Ts)),1);
%                 allIntratemp = unit(i).intraData(xrange,temp00)-ones(size(unit(i).intraData(xrange,temp00),1),1)*mean(unit(i).intraData(xrange,temp00),1)-10;
                allIntratemp = unit(i).intraData(xrange,temp00)-ones(size(unit(i).intraData(xrange,temp00),1),1)*mean(unit(i).intraData(1:int32(10e-3*Ts),temp00),1)+shiftmIntra0;
                                if bplotraw
                plot(xtime1, allIntratemp,'-k')
                                hold on
                                end                
                stemp0 = PSPParam(1,mIntra0,1/Ts,bIPSC,bIntra);    figure(figIND + k) %% get latency
                mIntra0 = mIntra0-temp0+shiftmIntra0;
                axis tight
                hold on
            end
            %% EX
            bIPSC =0;
            stemp1 = [0 0 0]; mIntra1 = 0;
            temp11 = intersect(tempIND,c1);
            if ~isempty(temp11)
                mIntra1 = mean(unit(i).intraData(:,temp11),2);
%                 sIntra1 =  std(unit(i).intraData(:,temp11),0,2);
                temp1 = mean(mIntra1(1:int32(10e-3*Ts)),1);
                allIntratemp = unit(i).intraData(xrange,temp11)-ones(size(unit(i).intraData(xrange,temp11),1),1)*mean(unit(i).intraData(1:int32(10e-3*Ts),temp11),1);
                if bplotraw
                plot(xtime1,allIntratemp,'-k')
                hold on 
                end
                stemp1 = PSPParam(1,mIntra1,1/Ts,bIPSC,bIntra);    figure(figIND + k)
                mIntra1 = mIntra1-temp1;
                plot(xtime1, mIntra1(xrange,:),'-r','LineWidth',2.5)
%                 plot(xtime1, mIntra1(xrange,:)-sIntra1(xrange,:),'--r','LineWidth',1.5)
 
                axis tight
                %                 ylim([unique(min(mIntra1))*3 unique(max(mIntra1))*2.5])
            end
            if ~isempty(temp00)
                plot(xtime1, mIntra0(xrange,:),'-c','LineWidth',2.5)
%                 plot(xtime1, mIntra0(xrange,:)-sIntra0(xrange,:),'--b','LineWidth',1.5)
            end
%             ylim([2.5*min([unique(min(mIntra0)) unique(min(mIntra1)) -10])  max([(unique(max(mIntra0))-shiftmIntra0)*2.5  40])+shiftmIntra0])


            %             ylim([-50+min(temp) 50+max(temp)]);
            if (subi==NN | ii == IMAX)
                ylabel('pA')
                xlabel('ms')
            else
                set(gca,'XTickLabel',' ')
            end
            
            %%% Calculate Latency from Mean
            title([stemp '    2PkEx:' num2str(stemp0(:,2),'%1.1f')  '    2PkInh:' num2str(stemp1(:,2),'%1.1f')],'Interpreter','none');
            set(gca,'Color','none','XGrid', 'on','XTick',[int32(min(xtime1):sGrid:max(xtime1))])
%             box off

            nSU =1+(subi-1)*NSU(1); % subplot in unit
%             sB = 1-((nSU-1)/(NS+NS*bH+eH*1.2)+(nSU-1)*bH +eH);
           sB = 1-((nSU)/(NS+NS*bH+eH*5)+(nSU-1)*bH );

            subplot('Position',[sL sB sW sH]) %% EE WAVEFOR<
            %% get waveform from electrode with largest spike
            mspike=mean(unit(i).spikes(:,tempIND),2)' ;
            ind_min = find(mspike==min(mspike),1,'first');
            electrode_min = floor(ind_min/(size(mspike,2)/size(EEgroup,2))) +1;
            %%
            temp1 = [-(ind_min-(size(mspike,2)/size(EEgroup,2))*(electrode_min-1))+1:abs((size(mspike,2)/size(EEgroup,2))*(electrode_min) - ind_min)];% index of min
            %set this at 0
            
            temp = [(electrode_min-1)*(size(mspike,2)/size(EEgroup,2))+1:(electrode_min)*(size(mspike,2)/size(EEgroup,2))];
            xtime = 1000*dt.*temp1;
            plot(xtime, unit(i).spikes(temp,tempIND),'-k')
            hold on
            plot(xtime,mspike(:,temp),'-r','LineWidth',2);
            axis tight
            xlim([min(xtime1) max(xtime1)])

             unit(i).SpikeParam = xSpikeParams2(mspike(:,temp),1/Ts*1000);
            spikeParam(i) = unit(i).SpikeParam(12) - unit(i).SpikeParam(11);
            title([sourcefile(1:end-4) '      ' 'A' sunitfilename ':' num2str(i) '    P2W:' num2str(spikeParam(i),'%1.1f') '   EX:' num2str(length(c1)) '   IN:' num2str(length(c0))],'Interpreter','none')
            set(gca,'Color','none','YTick',[],'XTickLabel',[' '],'XTick',[int32(min(xtime1):sGrid:max(xtime1))],'XGrid', 'on')
%             box off
        end
    end
    if breport
        dirtemp = 'REPORT';
        figdesc = [figString '_EE' num2str(triggerEE) '_W' num2str(k)];
        savefigure(writedirheader,dirtemp,figdesc,'emf')
    end

end

figure(570)
hist(spikeParam,20)
title('Histogram: Time From Trough to Peak')
xlabel('ms')

if breport
    figString = ['ASpikes_TtoP_Hist' sunitfilename];

    dirtemp = 'REPORT';
    figdesc = [figString '_EE' num2str(triggerEE)];
    savefigure(writedirheader,dirtemp,figdesc,'emf')
end
% 
% for k = 1:Nwindows
%     figure(figIND + k) %% Pair EE and Intra WINDOW
% 
%     if breport
%         dirtemp = 'REPORT';
%         figdesc = [figString '_EE' num2str(triggerEE) '_W' num2str(k)];
%         savefigure(writedirheader,dirtemp,figdesc,'emf')
%     end
% end
% if bsave
%         save(temp,'klusters','Center','num_c')  %% NOTE font size too big on axis (overlaps)
