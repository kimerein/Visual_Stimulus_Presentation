function [data dt time] = loadDAQData(LOADFILE,Chns,TrigRange)
% function [data dt time] = loadDAQData(LOADFILE,Chns,TrigRange)
% Load DAQ data using daqread and reshape into
% data = WAVEFORM x SWEEPS x SITES
% dt = sample interval
% time = time of each sample
%
% BA 091109

% NOTE: not considered adding sampRange suport, but data reshaping must be
% fixed to be compatible with this so it is commented out

bspecTrigger = 0;
% bspecSamp = 0;

if nargin>2    if ~isempty(TrigRange ); bspecTrigger = 1; end; end
% if nargin>3    if ~isempty(SampRange ); bspecSamp = 1; end; end

% if bspecTrigger && bspecSamp
% error ('The ''Samples'' and ''Triggers'' properties are mutually exclusive.')
% else
if      bspecTrigger
    [dataread time]  = daqread([LOADFILE '.daq'], 'Channels', Chns,'Triggers',TrigRange, 'Dataformat','native');
    % elseif      bspecSamp
    %     [dataread time]  = daqread([LOADFILE '.daq'], 'Channels', Chns,'Samples', SampRange, 'Dataformat','native');
else
    [dataread time]  = daqread([LOADFILE '.daq'], 'Channels', Chns, 'Dataformat','native');
end

dt = time(2) - time(1);
if nargout<3; clear time; end
dataread = single(dataread);

sweeplength = find(isnan(dataread),1,'first')-1; % length of each sweep
if isempty(sweeplength);
    sweeplength  = length(dataread); end
nsw = (length(dataread)+1)/(sweeplength+1); % number of sweeps/triggers
    
if ~nsw*(sweeplength+1)-1 == length(dataread)%check
    error('reshape dimensions wrong'); end
nchns = size(dataread,2);

% reshape to be WAVEFORM x SWEEPS x SITES
dataread = [dataread; nan(size(dataread,2),1)'];
if ~isinteger(nsw) % case where acquisition ended before the end of sweep
    nsw = floor(nsw);
    indEND = nsw*(sweeplength+1);    
else
    indEND = size(dataread,1);
end

data = reshape(dataread(1:indEND,:),sweeplength+1,nsw,nchns);
clear dataread;
data = data(1:sweeplength,:,:);


 
