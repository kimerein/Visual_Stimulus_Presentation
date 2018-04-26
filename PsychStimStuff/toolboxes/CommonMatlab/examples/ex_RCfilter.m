% ex_ RC filter
% from findspikes.m in % part of the Matteobox toolbox
cutoff = 1000;	% Hz. Not really. With 1000 it cuts off around 100 Hz.

% 1 - two RC circuits in series
deltat = 1/samplerate;
tau = 1/cutoff;
b = [ 1-(deltat/tau)^2 2*(deltat-tau)/tau ((deltat-tau)/tau)^2];
a = [ 1 2*(deltat-tau)/tau ((deltat-tau)/tau)^2];

% to look at the filter: freqz(b,a);

% -------------------- high-pass filter the potentials:
hipasdpots = zeros(nsamples,nstimuli);
% hipasdpots(:) = filter( b,a, potentials(:) ); 
% do it one by one. Filter is memory intensive.
for istim = 1:nstimuli
	hipasdpots(:,istim) = filter( b,a, potentials(:,istim) ); 
end
