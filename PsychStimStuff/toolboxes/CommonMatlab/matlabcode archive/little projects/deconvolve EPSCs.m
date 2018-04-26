N=12;
clear aa  %% use findOSC006.m to get Strough full of osc times
aa(1,:) = prodata(Strough(N)-10*1/output.dt/1000:Strough(N)+20*1/output.dt/1000)
LPFilter = 200;
dt =1e-4; % 10kHz
[B,A] = butter(2,2*LPFilter*dt,'low');
aaa =   filtfilt(B,A,aa);
aaa = aaa -mean(aaa(1:10));

figure(10)
plot(aaa);
aaa = [zeros(1,500) aaa zeros(1,500)]; %% pad so that the distribution can be long
%% aaa and EP need to globally defined. for use in Ferror
aaa = aaa(1:3:end);
LPFilter = 100;[B,A] = butter(2,2*LPFilter*dt,'low');
EP = filtfilt(B,A,uEPSC(:,2));
EP = EP(1:566);
EP = EP(1:3:end);
x0 = ones(1,abs(length(aaa)-length(EP))+1); %% first attempt defines length of distribution,D and is constrained to have certain length by EP and G
% x0 = [zeros(1,min(find(aaa>0))-1) ones; %% first attempt defines length of distribution,D and is constrained to have certain length by EP and G
A = diag(-1.*ones(length(x0),1));
b = zeros(1,length(x0));

option1 = optimset('MaxIter',1000,'Display','iter');
Result = fmincon(@Ferror,x0,A,b,[],[],0,50,[],option1);
figure(11)
clf
plot(Result*-100,'r')
hold on
plot(conv(Result,EP),'k')
plot(aaa,'b');



cc2 = conv(Result,EP);
plot(cc2)
hold on
plot(Result*-100,'r')
plot(aaa,'g')