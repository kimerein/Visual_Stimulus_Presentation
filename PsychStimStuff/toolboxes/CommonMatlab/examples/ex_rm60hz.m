%% example remove line frequency
% remove 60hz
params = 

          Fs: 10000
       fpass: [1 100]
      tapers: [3 5]
         pad: 2
         err: [1 0.0500]
    trialave: 1
          f0: []
%%

drm = rmlinesc(detrend(d(1:1/dt)),[],params,[],'y');