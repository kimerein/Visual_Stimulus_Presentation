%% Generate stimulus file
% W defines waveform to be put in
savefile = 'c:\DyClampStimulus_SynCond_10s.atf';

%% STIMULUS FILE PARAMETERS
MAXWAVEAMP = .9*10^4; % maximal amplitude of waveform
% should be set to match maximal allowed by output
% signal in clampex protocol
% e.g. 9000mV
%                     NOTE:(10V causes non-linear response either in
%                      DAC or CIA)
SwLength = 10000 ;   %duration in(ms) of generated stimulus file
dt = 1/11111 * 1000; % sampling interval (ms)

initOffset = 200; % offset from the begining of the stimulus file. after which waveforms can begin(ms)
t = [dt/1000:dt:SwLength]; % time in (ms)

Nj  = 100%% number of Sweeps to define
% Nj = N2;
%% CONDUCTANCE polarity i.e. biologicall this would be setting whether this
%% conductance is permeable to K or Na
% Direction of Currents
% Note Currents that are INTO the cell are NEGATIVE here:
% e.g. A conductance that reverses at 0mV will have dirn = -1;
dirn = -1;

%%%%%%%%%%%%%%%% GENERATE  WAVEFORM
%% this part creates
% a waveform that consists of synaptic waveform W repeated at an interval
% of intv.
% W is scaled randomly to 100% or 20% of MAXWAVEAMP

%% Synaptic W Definition inject
xW = [0:dt:25]; %% time course of waveform (ms)

% model of synaptic conduntance
tr = 2.5;
td =2.7;
W =    exp(xW/-td)-exp(xW/-tr) ; %see 1synapse dynamic clmap.cdr for description
W = W/max(W)*MAXWAVEAMP;
W = round(W)* dirn;


%% factors are used to scale Waveform to difference sizes
Lfactor = 1; % fraction of 10 that the large amplitude should be
Sfactor = 0.2;

%% W will be repeated at this interval
intv = 200; % interval to repeat W (ms)

%% declare final waveform variable that will be output into stimulus file
outWave = zeros(Nj,SwLength*1/dt,'double'); %% the output

nW = (length(outWave)-initOffset/dt)/(intv*1/dt); %% number of times to repeat W
while ((nW-1)*intv+length(W)*dt) > (SwLength-initOffset)%% make sure the last W is not cut off
    nW = nW-1;
end

n=0; m=0;
for j = 1:Nj
    bag = rand(1,int32(nW));
    for i = 1:int32(nW)
        %% randomly decide if large or small
        if bag(i) > 0.5
            amp = Sfactor*W;
            m = m+1;
        else
            amp = Lfactor*W;
            n = n+1;
        end
        outWave(j,1+((i-1)*intv+initOffset)*1/dt:((i-1)*intv+initOffset)*1/dt+length(W)) = amp;
    end
end
figure(1)
plot(t,outWave(1,:))
figure(2)
plot(outWave')

% n % number of large W
% m % number of small W

%%%%%%%%%%%%%%%%%%%%


%% out in Column vector form for output to stimulus file
outWave = double(outWave);
outdata = [t' outWave'];


%% SAVE to atf format
scol1 = 'Time (ms)';
scol2 = 'Sig (mV)';
fid = fopen(savefile,'w')
% fprintf(fid,'ATF 1.0\n0\t2\n"%s"\t"%s"\n',scol1,scol2); %%
fprintf(fid,'ATF 1.0\n0\t%d\n"%s"',Nj+1,scol1)
for j=1:Nj %% write colum titles
    fprintf(fid,'\t"%s"',scol2); %%
end
fprintf(fid,'\n');
for i=1:length(outdata)
    fprintf(fid,'%d',outdata(i,1)); %% must not save in %f or %g format so can't use save -ascii
    for j = 1:Nj
        fprintf(fid,'\t%d',outdata(i,j+1)); %% must not save in %f or %g format so can't use save -ascii
    end
    fprintf(fid,'\n'); %% must not save in %f or %g format so can't use save -ascii
end
fclose(fid)



