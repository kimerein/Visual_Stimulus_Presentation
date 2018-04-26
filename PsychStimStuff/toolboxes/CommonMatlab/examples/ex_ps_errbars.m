%% example PS with error bars 

params.err=[1 0.05]; [S,f,Serr]=mtspectrumc(detrend(df,'constant'),params);
% plot(f,10*log10(S),f,10*log10(Serr(1,:)),f,10*log10(Serr(2,:)));
plot(f,S);

%%
