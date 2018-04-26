
% bsubplotChn =1;% analyze LFP for each condition
% nChn = size(LFP,3);
% % nC = 4;
%%
aOFFSET = 0;

nC=length(C);
for i=1:nC; Cond(i).I_sw = find(TRIG==C(i)); end
%% plot individual sweeps
i = 1;
j = 1;
temp = Cond(C(i)).I_sw;
XL = [ 0 .5];
fid = 99; figure(99);clf;
subplot(5,1,1:4);
offsetplot([LFP(:,temp(1:3:end),j)],dt,3,fid);
xlim(XL);
subplot(5,1,5)
plotdata(mean(LFP(:,temp(1:end),j),2),dt/1000,'fid',fid);
xlim(XL)
%% plot average for each condition or Chn

% inSW_WINDOW = [.8 1.5]/dt;
% average LFP
xtime = [1:size(LFP,1)]*dt;
fid = 21;    figure(fid);clf; clear sleg lh h mLFP;
for i = 1:nC % for each condition
    mLFP = [];
    for j = 1 :nChn % for each chn
        %         (i-1)*size(LFP,3)+j
        
        figure(fid);
        % subplot each Chn
        if bsubplotChn
            colorI = i;
            h(j) = subplot(nChn,1,j);          cmap = jetm(nC);
            sleg{i} = num2str(C(i));
        else % subplot each cond
            h(i) = subplot(nC,1,i);            cmap = jetm(size(LFP,3));
            colorI = j;
            sleg{j} = num2str(j);
        end
        
        if 0      %plot all traces
            %             lh(colorI)  = mplot(xtime, squeeze(LFP(:,Cond(i).I_sw,j))', 'Color', cmap(colorI,:));
            lh(colorI)  = plot(xtime, squeeze(LFP(:,Cond(C(i)).I_sw,j))');
            set(lh(colorI) , 'ButtonDownFcn', {@raise_me, lh(colorI)});hold on;
            
            
        else % just plot mean
            % subplot each Cond
            mLFP(j,:) = mean(LFP(:,Cond(i).I_sw,j),2);
            %            lh(colorI) = plot(xtime,mLFP(j,:));
            lh(colorI) = plot(xtime,mLFP(j,:)+aOFFSET*i,'color',cmap(colorI,:),'linewidth',2); hold on;
            set(lh(colorI) , 'ButtonDownFcn', {@raise_me, lh(colorI)});hold on;
            
            %             plot(xtime,std(LFP(:,Cond(C(i)).I_sw,j),1,2),'color',cmap(colorI,:),'linewidth',1,'linestyle','--'); hold on;
            %                          plot mean and std with area
            %                                       [lh,ph] = errorarea(xtime,mean(LFP(:,Cond(C(i)).I_sw,j),2), std(LFP(:,Cond(C(i)).I_sw,j),1,2)); hold on;
            %                             set(lh, 'Color', brighten(cmap(colorI,:), -0.6), 'ZData', repmat(j, size(get(lh,'XData'))));
            %                             set(ph, 'FaceColor', cmap(colorI,:), 'ZData', repmat(j, size(get(ph,'XData'))), 'FaceAlpha', 0.8);
            
        end
        axis tight;  if j == nChn;xlabel('time (s)');end
        % plotset(1);
    end
    
    if bCSD % compute CSD
        figure(fid+5);
        computeCSD(mLFP,1/dt)
    end
end
% linkprop(h,{'box','TickDir','Xlim'});
legend(lh,sleg)
%% trial by trial coherence between sweeps

i = 1;% condition
j = 1 ;% chn
scolor = 'k';
bphase = 0;
bpw = 0;
cohparam.Fs = 1/dt;
cohparam.fpass = [0 200]; %kHz
cohparam.err = [0 0.05];
cohparam.tapers = [3 5];
% bplotlog = 0;
% win = size(LFP,1);
%
% for i = 1:size(LFP,2)
%      ind = find([1:size(LFP,2)]~=i);
%      for k = 1:length(ind)
% %         [CxyCohere(:,k),phi(:,k),S12(:,k),S1(:,k),S2(:,k),F]=coherencysegc(squeeze(LFP(:,ind(k),j)),squeeze(LFP(:,i,j)),win,cohparam);
%      end
%      cxy(:,i) = mean(CxyCohere,2);
% end
% TWIN = round([0.2 0.4]/dt);
% TWIN2 = round([1.2 1.4]/dt);
TWIN = round([1.0 1.2]/dt);
TWIN2 = round([.6 .8]/dt);
% get an estimate by looking for coherence between two different sweeps look at mean
clear  Cxy Cxy2 phi S12 S1 S2;

for j = 1: size(LFP,3) % average across channels
    %        [CxyCohere,phi,S12,S1,S2,F]=coherencyc(squeeze(LFP(1:end/2,1:end-1,j)),squeeze(LFP(1:end/2,2:end,j)),cohparam);
    [tempC phi,S12,S1,S2,F]=coherencyc(squeeze(LFP(TWIN(1):TWIN(2),1:end,j)),squeeze(LFP(TWIN(1):TWIN(2),end:-1:1,j)),cohparam);
    Cxy(j,:) = mean([tempC],2);
    [tempC,phi,S12,S1,S2,F]=coherencyc(squeeze(LFP(TWIN2(1):TWIN2(2),1:end,j)),squeeze(LFP(TWIN2(1):TWIN2(2),end:-1:1,j)),cohparam);
    Cxy2(j,:) = mean([tempC],2);
end

if j > 1
    Cxy = mean(Cxy,1);
    Cxy2 = mean(Cxy2,1);
end
figure(98);clf;
plot(F,Cxy,'-r'); hold all;
plot(F,Cxy2,'-k');

title('coherence between trials')
%% Power spectrum
sparams.Fs = 1/dt; sparams.fpass = [0 100]; sparams.err = [2 0.05]; sparams.trialave = 0;
pparam.logplot = 1;
fid = 50; figure(fid);clf;
% XWIN = [round(.8/dt):round(1/dt)];
XWIN = [TWIN(1):TWIN(2)];
% plotps(mean(LFP(XWIN,:),2)',dt,fid,[0 100],sparams);
pparam.color = [1 0 0];
plotpsCHX(mean(LFP(XWIN,:),2),dt,fid,[0 100],sparams,pparam);
hold all;
% XWIN = [round(1.2/dt):round(1.4/dt)];
XWIN = [TWIN2(1):TWIN2(2)];
% plotps(mean(LFP(XWIN,:),2)',dt,fid,[0 100],sparams);
pparam.color = [0 0 0];
plotpsCHX(mean(LFP(XWIN,:),2),dt,fid,[0 100],sparams,pparam);

%%
% % SPectrogram
% fid = 30;    figure(fid);clf;
% if 1
%     Res = 2; nfft = 2^nextpow2((1/dt/2)/Res);   % define resolution of fft
%     clear sparams
%     movingwin = [0.2 0.03]; %window, windowstep
%     sparams.Fs = 1/dt; sparams.fpass = [0 100]; sparams.err = [2 0.05]; sparams.trialave = 0;
%
%
%     for i = 1:nC % for each condition
%         for j = 1:size(LFP,3) % for each chn
%             %         (i-1)*size(LFP,3)+j
%             subplot(size(LFP,3),nC,(i-1)*size(LFP,3)+j)
%             % spectrogram
%             [S,t,f,Serr]=mtspecgramc(mean(LFP(:,Cond(C(i)).I_sw,j),2)',movingwin,sparams);
%             plot_matrix(S,t,f,'n')
%         end
%     end
%
% end
%
% % for all conditions
% for j = 1:size(LFP,3) % for each chn
%
%     subplot(size(LFP,3),1,j)
%     % spectrogram
%     [S,t,f,Serr]=mtspecgramc(mean(LFP(:,:,j),2)',movingwin,sparams);
%     plot_matrix(S,t,f,'l')
% end
%
%
% %% Coherence
% ifid = 30;    figure(fid);clf;
% if 1
%     Res = 2; nfft = 2^nextpow2((1/dt/2)/Res);   % define resolution of fft
%     clear sparams
%     movingwin = [0.05 0.01]; %window, windowstep
%     sparams.Fs = 1/dt; sparams.fpass = [0 200]; sparams.err = [2 0.05]; sparams.trialave = 0;
%
%     %      coherence of one LFP with all others
%     for i = 1:nC % for each condition
%         for j = 2:size(LFP,3) % for each chn
% %             subplot(size(LFP,3)-1,nC,(i-1)*(size(LFP,3)-1)+j)
%             % cohere
%             %             plotcohgramCHX(df,df2,dt,movingwin,sparams,blogplot)
%             plotcohsegCHX(squeeze(LFP(:,Cond(C(i)).I_sw,1))',squeeze(LFP(:,Cond(C(i)).I_sw,j))',dt,[1 size(LFP,1)],'fid',ifid,'color', cmap(j,:));hold on;
%             plotcohsegCHX(squeeze(LFP(round(.1/dt):round(.3/dt),Cond(C(i)).I_sw,1))',squeeze(LFP(:,Cond(C(i)).I_sw,j))',dt,[1 size(LFP,1)],'fid',ifid,'color', cmap(j,:));hold on;
%         end
%     end
%
% end
%
% % clear sleg
% % ifid = 31;
% % XWIN = [round(.1/dt):round(.3/dt)];
% %     for j = 2:size(LFP,3) % for each chn
% % %             subplot(size(LFP,3)-1,nC,(i-1)*(size(LFP,3)-1)+j)
% %             % cohere
% %             %             plotcohgramCHX(df,df2,dt,movingwin,sparams,blogplot)
% %             plotcohsegCHX(squeeze(LFP(XWIN,Cond(C(i)).I_sw,1))',squeeze(LFP(XWIN,Cond(C(i)).I_sw,j))',dt,[length(XWIN)]*dt,'fid',ifid,'color', cmap(j,:));hold on;
% %     end
% %
% %     ifid = 32;
% % XWIN = [round(1.4/dt):round(1.6/dt)];
% %         for j = 2:size(LFP,3) % for each chn
% % %             subplot(size(LFP,3)-1,nC,(i-1)*(size(LFP,3)-1)+j)
% %             % cohere
% %             %             plotcohgramCHX(df,df2,dt,movingwin,sparams,blogplot)
% %             plotcohsegCHX(squeeze(LFP(XWIN,Cond(C(i)).I_sw,1))',squeeze(LFP(XWIN,Cond(C(i)).I_sw,j))',dt,[length(XWIN)]*dt,'fid',ifid,'color', cmap(j,:));hold on;
% %             sleg{j} = num2str(j);
% %         end
% %         sleg{1} = ''
% %         legend(sleg)
% %
% %         ifid = 41;figure(ifid);clf;
% % XWIN = [round(.1/dt):round(.3/dt)];
% %     for j = 2:size(LFP,3) % for each chn
% % %             subplot(size(LFP,3)-1,nC,(i-1)*(size(LFP,3)-1)+j)
% %             % cohere
% %             %             plotcohgramCHX(df,df2,dt,movingwin,sparams,blogplot)
% %             plotcohsegCHX(squeeze(LFP(XWIN,Cond(C(i)).I_sw,5))',squeeze(LFP(XWIN,Cond(C(i)).I_sw,j))',dt,[length(XWIN)]*dt,'fid',ifid,'color', cmap(j,:));hold on;
% %     end
% %
% %     ifid = 42;figure(ifid);clf;
% % XWIN = [round(1.5/dt):round(1.7/dt)];
% %         for j = 2:size(LFP,3) % for each chn
% % %             subplot(size(LFP,3)-1,nC,(i-1)*(size(LFP,3)-1)+j)
% %             % cohere
% %             %             plotcohgramCHX(df,df2,dt,movingwin,sparams,blogplot)
% %             plotcohsegCHX(squeeze(LFP(XWIN,Cond(C(i)).I_sw,5))',squeeze(LFP(XWIN,Cond(C(i)).I_sw,j))',dt,[length(XWIN)]*dt,'fid',ifid,'color', cmap(j,:));hold on;
% %             sleg{j} = num2str(j);
% %         end
% %         sleg{1} = ''
% %         legend(sleg)
%
% % % for all conditions
% % for j = 1:size(LFP,3) % for each chn
% %
% %     subplot(size(LFP,3),1,j)
% %     % spectrogram
% %     [S,t,f,Serr]=mtspecgramc(mean(LFP(:,:,j),2)',movingwin,sparams);
% %     plot_matrix(S,t,f,'l')
% % end
