function [STA Nspikes Param]= getSTA(FilenameSpkSort,UNIT,DAQchn,STF,STAparams)
% function [STA MovieName]= getSTA(FilenameSpkSort,UNIT,DAQchn,STF,STAparams)
%
% INPUT
%     - FilenameSpkSort
%     - UNIT  scalar indicating unit to analyze
%     - DAQchn to use to extract frame times
%     - STF (optional) struct with field "filename" containing DAQ to analyze
%             % if not supplied will use all DAQfile in FilenameSpkSort that
%             have Vstim movie = 1st DAQfile with Vstim movie)
%     - STAparams (optional) struct 
%            .tbin  compute STA in this size step (ms) <5>
%            .WOI  [before_spk after_spk] compute STA over this interval (ms) <[200 50]> , i.e. 200 ms before spike to 50 ms after
%            spike
%
% BA 011910

DAQSAVEPATH = []; SORTEDSPIKES_SAVEPATH = []; LOGFILEPATH = [];
MOVIE_PATH = [];
s_RIGSPECIFIC_SPIKESORT_CONFIG
     
if exist('STAparams','var')
    tbin = STAparams.tbin;
    WOI = STAparams.WOI;
else % defaults
    tbin = 5/1e3; % ms
    WOI = [200 50]/1e3; % ms
end

load(fullfile(SORTEDSPIKES_SAVEPATH,FilenameSpkSort),'spikesortdata','rawdatainfo');
dt = rawdatainfo.dt;
% analyze only STF.filenames if specified otherwise try to use all
% rawdatainfo.filenames
if exist('STF','var') & ~isempty(STF); analyzefilenames = {STF.filename}; else analyzefilenames = {rawdatainfo.filename}; end 

%% get stimulus in form Time x Space
MovieName = [];
maskismovie = ones(size(analyzefilenames));
for i = 1:length(analyzefilenames)
    [junk Param] = getPsychStimParameters(analyzefilenames{i},DAQSAVEPATH, LOGFILEPATH);
    if isempty(findstr(Param.StimulusName, 'Movies'))
        maskismovie(i) = 0;
        if exist('STF','var') % if filename was specified by user and not movie throw error
            s = sprintf('PyschStimParameters associated with %s do not contain STA movie stimuluserror',analyzefilenames{i}); error(s);
        end
    else
        if isempty(MovieName)
            MovieName = Param.MovieName;
        elseif maskismovie(i) == 0; ~isequal(MovieName,Param.MovieName); printf('%s has a different movie: %s',analyzefilenames{i}, Param.MovieName); end;        
    end
end


analyzefilenames = analyzefilenames(find(maskismovie));
 
% this takes a couple secs probably not worth saving
load(fullfile(MOVIE_PATH,MovieName))
sz = size(moviedata);
Stim = reshape(moviedata,sz(1)*sz(2),sz(3));
EFRAMES = size(moviedata,3); % expected number of frames
clear moviedata;


%% get frame times in struct with each sweep
% reformat for function 

tempSTF = cell2struct(analyzefilenames, 'filename', length(analyzefilenames));
timeOfFrame = getDAQFramesTimes(DAQSAVEPATH,tempSTF,EFRAMES,DAQchn);

% spiketimes into convert into struct with each sweep

fileind = find(ismember({rawdatainfo.filename},{tempSTF.filename}));
totalspikes = 0;
for i = 1:length(analyzefilenames);
    fileind = find(ismember({rawdatainfo.filename},tempSTF(i).filename)); % find files in sort that match files that need to be analyzed (i.e. in tempSTF) and keep in tempSTF order
    if ~isempty(fileind)
        if length(unique(spikesortdata.file(fileind).sweep))>1; warning('DAQfile appears to have multiple sweeps.\n ONLY First sweep used'); end % 1 per file, can be fixed later
        mask = (spikesortdata.file(fileind).sweep  == 1)& (spikesortdata.file(fileind).assign  == UNIT);
        swdata(i).spiketimes = spikesortdata.file(fileind).spiketimes(mask)*dt;
        totalspikes = totalspikes + length(swdata(i).spiketimes);
    else
        s = sprintf('%s does not exist in spksortFN: %s ', tempSTF(i).filename, FilenameSpkSort); error(s); end  
end


WOIbin = WOI/tbin; % convert to bins
ST_framesNum = zeros(totalspikes,sum(WOIbin)+1,'int16');
for swN = 1:length(swdata)
    lstTime = timeOfFrame(swN).frametime(end) + WOI(1); % time after which not to accept spikes (because spikes may occur after frames stop being shown)
    
    % convert frametimes to bins (super sampling) TODO vectorize (not actually slow)
    frame_bin = floor((timeOfFrame(swN).frametime)/tbin)+1; % tbin that each frame is in
    % there is no zero tbin
    resample_frame_bin = zeros(round(lstTime/tbin),1,'int32'); %sweep length vector one element for each bin
    for frameN = 1:length(frame_bin)-1 % for each frame
        resample_frame_bin(int32([frame_bin(frameN):(frame_bin(frameN+1)-1)]))= frameN; %
    end
    resample_frame_bin(frame_bin(end):end) = length(frame_bin);
    
    
    stime_bin = round(swdata(swN).spiketimes/tbin);%     round stimes to WOIstep bins
    % throw away spike times without window before or after
    stime_bin = stime_bin( (stime_bin>WOIbin(1)) & ((stime_bin + WOIbin(1))<lstTime/tbin));
    [insertrow temp] = find(ST_framesNum==0,1,'first');
    temp = getWOI(resample_frame_bin,stime_bin,WOIbin); %
    ST_framesNum(insertrow:insertrow+size(temp,1)-1,:)= temp; % get chunks of frames triggered on spike
end
ST_framesNum = ST_framesNum';
clear I J stime_bin frame_bin lstTime resample_frame_bin temp timeOfFrame 
%% convert ST_Frames into real stimlus
temp = ST_framesNum>0; % get rid of zerosi
ST_frames= Stim(:,ST_framesNum(temp));
Nspikes = sum(temp(1,:)>0);
ST_frames = reshape(ST_frames,size(ST_framesNum,1)*size(Stim,1),Nspikes);

clear ST_framesNum;

STA = mean(ST_frames,2);
STA = reshape(STA,sz(1)*sz(2),(sum(WOIbin)+1));

if 0 % whiten
    CovStim = cov(single(Stim'));
    InvCovStim = inv(CovStim);
    
    for i = 1:size(STA,2)
    wSTA(:,i) = InvCovStim*STA(:,i);
    end
    wSTA = reshape(wSTA,sz(1),sz(2),(sum(WOIbin)+1));

end
clear Stim;

STA = reshape(STA,sz(1),sz(2),(sum(WOIbin)+1));
% whitenning doesn't leave anything, maybe beceause of warning from inverting
% covariance matrix?
% PLAY STA BACK

%% STC (from Rust tutorial)
% % take out subset of STA and ST_Frames
% selectbin = 20; % select time bin desired
% subSTA = single(STA(:,:,selectbin));
% subST_frames = single(ST_frames((selectbin-1)*sz(1)*sz(2)+1:selectbin*sz(1)*sz(2),:));
% 
% clear STA ST_frames
% 
% ST_frames2 = subST_frames'-subST_frames'*subSTA(:)*subSTA(:)';
% % NOT SURE is right
% 
% covv = ST_frames2'*ST_frames2./(size(subST_frames,2)-1);
% 
% subplot(2,2,1);
% % The diagonal of this matrix describes the variance, which swamps the
% % covariance between dimensions. To visualize the covariance matrix, we
% % will set the diagonal to zero so we can see the covariance structure:
% imagesc(covv.*(1-eye(size(covv,1))));colormap(gray);
% title('Covariance matrix')
% 
% [evect,eval] = sortedEig(covv);  eval = diag(eval); %% takes forever
% 
% subplot(2,2,2)
% mxx = max(max(eval(1:end-1))).*1.3;
% plot(eval(1:end-1),'.')
% axis([0 length(evect) 0 mxx])
% xlabel('rank number')
% ylabel('eigenvalue (variance)')
% title('PCA of the covariance matrix')
% 
% % plot filters with large or small variance
% % is the first evect the STA? shouldn't it be projected out?)
%  figure(2);
%  imagesc(reshape (evect(:,1),sz(1),sz(2)))
%  
% % whiten
% 
% % covariance matrix of stimulus
% 

