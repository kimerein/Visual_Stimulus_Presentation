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
    xf = [0:0.1:size(x)]';
    f = polyval(p,xf);
end
figure(3)
plot(V(:,:),I(:,:))

figure(4)
% plot([-1:0.1:0.9],p(:,1))
plot(p(:,1))
title('Max 2nS , 400pA/V')
ylabel('Conductance (nS)')


