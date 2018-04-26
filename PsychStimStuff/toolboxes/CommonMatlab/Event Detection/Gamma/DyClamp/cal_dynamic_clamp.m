% output = importStruct_abf('C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\dyclamp_test.abf',-1,0); 
%% 2nA/V in this file

% output = importStruct_abf('C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\dyclamp_test_400pA_V.abf',-1,0);
output = importStruct_abf('C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\2006_10_22_0016.abf',-1,0);
I = output.data(160:4000,3:3:end);
V = output.data(160:4000,2:3:end);
%% inport file containing steps in conductance
for i=1:size(V,2)
    x = V(:,i);
    y = I(:,i);
    p(i,:) = polyfit(x,y,1);
end
figure(3)
plot(V(:,:),I(:,:))

figure(4)
xp = [-9:1:9]; %% assumes 10nS max
 plot(xp,p(:,1),'.b')
cal = polyfit(xp',p(:,1),1); %% find calibration factor
%% NOTE it is best to manually remove low conductance values when making
%% this fit because they throw off the fit
hold on
    xf = [1.1*min(xp):0.1:max(xp)*1.1]';
    f = polyval(cal,xf);
    plot(xf,f,'-k')
title(['Max 10nS , (400pA/V)' ' Real Cond = ' num2str(cal(1),'%1.2f') '* CIA Cond + '  num2str(cal(2),'%1.2f')]  )
ylabel('Conductance (nS)')



