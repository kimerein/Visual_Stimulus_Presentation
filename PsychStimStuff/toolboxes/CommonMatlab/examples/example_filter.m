
[B,A] = cheby2(9,20,.5,'high');
HIGHCUT = 100;%Hz
[B,A] = butter(2,2*HIGHCUT/Ts,'high'); 

figure(3)
freqz(B,A,512,1000) % plot filter
% zplane(Hd)
% f = dfilt.filter(b,a,single(output.data(:,6+105*N_chn)))
Hd = dfilt.df2t(B,A);
% f = filter(Hd,single(output.data(:,6+105*N_chn)))

f = single(filter(Hd,single(allElectrodesXspikes(:,:))))

% close 1
figure(1);
plot(f)
figure(2)
% plot(output.data(:,6+105*N_chn)+5000,'-r')
plot(allElectrodesXspikes(:,1),'-r')

N = 2;
Y = filtfilt(B,A,double(a));
figure;
plot(a(:,N)-mean(a(:,N)),'r')
hold on
plot(Y(:,N))
