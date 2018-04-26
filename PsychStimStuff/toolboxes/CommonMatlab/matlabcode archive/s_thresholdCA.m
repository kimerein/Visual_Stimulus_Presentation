%% analyze ROIs extracted from imagej
% pasted by hand Z data
% 

%%
% plot(a(:,2:end))
a2 = a - repmat(mean(a),size(a,1),1);
a2 = a2(:,2:end);
% plot(a2')

dt = 80e-3; %sec
%% threshold Ca transients
sa = mean(std(a2));

th = sa*1.2;
indx = threshdiff(a2(1:end)*-1,th*-1);
% indx = thresh(a2(1:end)*-1,th*-1);

event = getWOI(a2,indx,[-1 8]);
%% plot rois and thresholdcrossing

fid = 31;
figure(fid)
clf;
off = offsetplot(a2,1,5,fid);
hold on
rasterplot(indx,1,length(a2),off,fid)

% figure(2)
% plot(event(:,:)')
% 
figure(5)
plot(a2(:,3))
hold all
plot(a2(:,18))

%%
swl = size(a2,2);
figure(10)
clf
plot(a2(1:end))
hold on
plot((a2(1:end)*-1<th*-1)+15,'-r')
plot(indx,ones(size(indx))*15,'xg')
%% are the event transients like qi's?


%% OTHER
%% high pass filter
for i=1:size(a2,2)
    af(:,i) = filterdata(a2(:,i),dt,1/dt/3,0);
    af(:,i) = filterdata(af(:,i),dt,1/dt/100,1);
end
figure(21)
clf
plot(a2(1:end))
hold all
plot(af(1:end)+50)


%% look for bleaching
for i=1:size(a2,2)
    af(:,i) = filterdata(a2(:,i),dt,1/dt/100,0);
%     af(:,i) = filterdata(af(:,i),dt,1/dt/100,1);
end
figure(22)
clf
plot(a2(1:end))
hold all
plot(af(1:end)+50)

figure(20)
clf
plotdata(a2,dt,'chn',[1 5 10 12],'fid',20)
hold all
plotdata(af,dt,'chn',[1 5 10 12],'fid',20)