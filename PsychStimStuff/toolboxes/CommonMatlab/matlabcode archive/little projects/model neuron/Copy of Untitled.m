    scaleI=10;
    scaleU=25;
    cells=neighbors(1); %total-2 ;%extneigh(1);
    time=1:ceil(((i-1)/tausPerRecSamp));

% figure;    plot([lts(cells,time);ltsI(cells,time)/scaleI; ltsU(cells,time)/scaleU ]');
figure(100);     
% plot(time.*tau,lts(extneigh(1),1:size(time,2)),'r');hold on;
plot(time.*tau,lts(cells,1:size(time,2)),'b');title(stemp);
% 
% % plot(time.*tau,lts(:,1:size(time,2)),'b');title(stemp);
% figure(102);     
% plot(time.*tau,lts(extneigh(1),1:size(time,2)),'r');hold on;title(stemp);
% % plot(time.*tau,lts(cells,1:3331),'b');
% 
% cells = unique(firings(:,2));
% figure(103);
% for(i=1:size(cells))
%     plot(time.*tau,lts(cells(i),1:size(time,2))-50*i,'b');title(stemp);
%     hold on;
% end
% % spike = diff(lts'>0)'>0; % find spikes
% % % spike = lts'>20; % find spikes
% % % A = xcorr(spike(neighbors(2),:)',spike(neighbors(2),:)',length(spike)/2);
% % 
% % % doesn't work why?
% % figure(101)
% % title(stemp)
% % B = xcov(spike(neighbors(2),:),spike(neighbors(3),:));
% % plot(B);
% % hold on
% % ii=1;
% % % for ii=1:length(extneigh)
% %     C = xcov(spike(neighbors(2),:),spike(extneigh(ii),:));
% %     plot(C -5 ,'r');
% % %     hold onun
% % % end
% % 
% % % FFT why no clear peak 
% % % y1 = fft(spike(neighbors(2),1:end/2));
% % % m1 = abs(y1);  % Magnitude 
% % % f1 = (0:length(y1)-1)*99/length(y1); % Frequency vector
% % % %% FFT why zero
% % % y2 = fft(spike(neighbors(2),end/2:end));
% % % m2 = abs(y2);  % Magnitude 
% % % f2 = (0:length(y2)-1)*99/length(y2); % Frequency vector
% % % figure(99);plot(f1,m1);hold on; plot(f2,m2,'r'); title('Magnitude');
% % 
% % % a =mean(spike(extneigh(:),:));
