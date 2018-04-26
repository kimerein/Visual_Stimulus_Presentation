% example powerspectrum

% useing mtm

    maxres = 1;
    nw = maxres*(length(x)*dt)
    [pyy,Pxxc,f] = pmtm(double(x),{nw,'trace'},[],1/dt);
    hpsd = dspdata.psd([Pxx Pxxc],'Fs',fs);
plot(hpsd)
