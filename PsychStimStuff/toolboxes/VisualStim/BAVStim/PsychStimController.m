function varargout = PsychStimController(varargin)
% PSYCHSTIMCONTROLLER M-file for PsychStimController.fig
%      PSYCHSTIMCONTROLLER, by itself, creates a new PSYCHSTIMCONTROLLER or raises the existing
%      singleton*.
%
%      H = PSYCHSTIMCONTROLLER returns the handle to a new PSYCHSTIMCONTROLLER or the handle to
%      the existing singleton*.
%
%      PSYCHSTIMCONTROLLER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PSYCHSTIMCONTROLLER.M with the given input arguments.
%
%      PSYCHSTIMCONTROLLER('Property','Value',...) creates a new PSYCHSTIMCONTROLLER or raises the
%      existing singleton*.  Starting from the left, property value pairs
%      are
%      applied to the GUI before PsychStimController_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PsychStimController_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help PsychStimController

% Last Modified by GUIDE v2.5 19-Aug-2009 18:52:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @PsychStimController_OpeningFcn, ...
    'gui_OutputFcn',  @PsychStimController_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before PsychStimController is made visible.
function PsychStimController_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PsychStimController (see VARARGIN)

rigSpecific;         %%% load file with rig-specific parameters
%% TO DO add case where rigSpecific doesn't exist

handles.USER.moviedirpath  = PSC_moviedirpath;
handles.USER.paramdirpath= PSC_paramdirpath;
handles.USER.logdirpath= PSC_logdirpath;
handles.USER.DAQPCIP = PSC_DAQ_PC_IP;
handles.USER.GAMMATABLE = PSC_GAMMATABLE;

handles.USER.Sync.kNone = 1;
handles.USER.Sync.ktdtSync = 2;

% Choose default command line output for PsychStimController
handles.output = hObject;
screennums = Screen('Screens'); % if there is just one monitor use that unstaed of th default monitor
if length(screennums)==1; % checkk if more than 1 monitor exists
    set(handles.ScreenNum,'String',num2str(screennums))
end
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PsychStimController wait for user response (see UIRESUME)
% uiwait(handles.figure1);

ScreenNum_Callback(handles.ScreenNum,eventdata,handles);
StimType_Callback(handles.StimType,eventdata,handles);
Var1_Callback(handles.Var1,eventdata,handles);
Var2_Callback(handles.Var2,eventdata,handles);


% --- Outputs from this function are returned to the command line.
function varargout = PsychStimController_OutputFcn(hObject, eventdata, handles) %#ok
varargout{1} = handles.output;


% --- Executes on selection change in StimType.
function StimType_Callback(hObject, eventdata, handles)

StimType = get(hObject,'Value');
%%% disable fields that aren't appropriate to this stimulus type

if StimType == 1 || StimType == 6     %% drifting or counterphase gratings
    set(handles.Duration,'Enable','on');
    set(handles.Speed0,'Enable','off');
    set(handles.TempFreq0,'Enable','on');
    set(handles.Phase0,'Enable','off');
end

if StimType == 2     %% drifting bars
    ScreenSizeDegX = str2double(get(handles.SizeX,'String')) * ...
        atan(1/str2double(get(handles.ScreenDist,'String'))) * 180/pi;
    Duration = ScreenSizeDegX/str2double(get(handles.Speed0,'String'));
    set(handles.Duration,'String',num2str(Duration));
    set(handles.Duration,'Enable','off');
    set(handles.Speed0,'Enable','on');
    set(handles.TempFreq0,'Enable','off');
end

if StimType==6 %% counterphase gratings
    set(handles.Phase0,'Enable','on');
end

if StimType == 1 || StimType == 2 || StimType == 6 %% drifting bars, drifting or counterphase gratings
    set(handles.Orient0,'Enable','on');
    set(handles.Freq0,'Enable','on');
    set(handles.Contrast0,'Enable','on');
    set(handles.SelectMovieName,'Enable','off');
    set(handles.MovieName,'Enable','off');
    set(handles.MovieMag,'Enable','off');
    set(handles.MovieRate,'Enable','off');
    set(handles.phasePeriod,'Enable','off');
    set(handles.stimulusGroups,'Enable','off');
end

if StimType == 7 %% spot
    set(handles.PositionX0,'Enable','on');
    set(handles.PositionY0,'Enable','on');
    set(handles.Duration,'Enable','on');
end

if StimType == 3 %% movie
    set(handles.Orient0,'Enable','off');
    set(handles.Speed0,'Enable','off');
    set(handles.Freq0,'Enable','off');
    set(handles.Contrast0,'Enable','off');
    set(handles.PositionX0,'Enable','on');
    set(handles.PositionY0,'Enable','off');
    set(handles.Length0,'Enable','off');
    set(handles.Duration,'Enable','on');
    set(handles.SelectMovieName,'Enable','on');
    set(handles.MovieName,'Enable','on');
    set(handles.MovieMag,'Enable','on');
    set(handles.MovieRate,'Enable','on');
    set(handles.phasePeriod,'Enable','on');
    set(handles.stimulusGroups,'Enable','on');
    set(handles.TempFreq0,'Enable','off');
    set(handles.Phase0,'Enable','off');
end

if StimType == 3
    % set wait interval to default 0 -- it's too easy to screw this up when
    % showing movies, causing unusual contimage phase behavior -- this is a
    % hack but I don't have better ideas MSC
    set(handles.WaitInterval,'String','0');
end

Var1_Callback(handles.Var1,eventdata,handles);
Var2_Callback(handles.Var2,eventdata,handles);


% --- Executes during object creation, after setting all properties.
function StimType_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to StimType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on button press in RunBtn.
function RunBtn_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to RunBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try %% put this is a try/catch, so that any crash won't leave screen hung


    clear mex
    %%% save parameters automatically
    cd(handles.USER.logdirpath);
    if ~exist(date,'dir')
        s = sprintf('mkdir %s',date);
        dos(s);
    end

    fname = fullfile(date,datestr(clock,30));
    SaveParams(handles,fname);

    %%% BA send stimuls filename to other computer
    if ~isempty(handles.USER.DAQPCIP)
    u = udp(handles.USER.DAQPCIP,9093,'LocalPort',9094);
    fopen(u);
    end

    InitializeMatlabOpenGL;   %%%necessary for OpenGL calls (like ClutBlit)

    %%% display description
    Duration = str2double(get(handles.Duration,'String'));
    FrameHz = round(str2double(get(handles.FrameHz,'String')));
    whichScreen = str2double(get(handles.ScreenNum,'String'));
    [window,windowRect]=Screen(whichScreen,'OpenWindow',0);   %%% open grey window
    Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); % BA needed for masking (could find another way )


    white = WhiteIndex(window);
    black = BlackIndex(window);
    grey = round(0.5*(black+white));

    imageRect = windowRect;
    Screen('DrawText',window,sprintf('Generating stimuli'),10,30,white);
    Screen('Flip',window);

    ScreenSizeDegX = str2double(get(handles.SizeX,'String'))*atan(1/str2double(get(handles.ScreenDist,'String')))*180/pi;
    degPerPix = ScreenSizeDegX/windowRect(3);

    handles.degPerPix = degPerPix;
    params.nCond = size(handles.orient,2);
    stim = get(handles.StimType,'Value');

    % nReps is the number of repetitions of the entire stimulus sequence
    % fractional values of nReps are useful to truncate movies
    nReps = str2double(get(handles.nReps,'String'));
    if nReps <= 0 %% loop almost-infinitely
        nReps = 10000000;
    end

    nStimulusRepetitions = 1; % defaults to 1 for non-movie stimuli

    %%%for clut animation
    if stim == 1 || stim == 2 || stim == 4 || stim == 5 || stim == 6 || stim == 7
        clut=1;
        sizeLut=256;
        offclut=zeros(sizeLut,3);
        offclut(:,:)=grey;   %default
    else %% stim == 3
        clut=0; % movie
    end

    switch stim
        %% TO DO BA added box in corner of screen that flips from gray to white
        %% whenever fram changes
        %%%%%%  drifting and counterphase gratings %%%%%
        case {1,6}
            % screen_gamma=2;
            textures = zeros(params.nCond,1);
            for c = 1:params.nCond
                if get(handles.StimType,'Value')==1   %%% drift vs counterphase
                    condNAME = 'Drifting Gratings';% BA

                    %% NOTE %%% THIS clut METHOD introduces weird artifacts on the
                    %% screen,(should remove it from other stimuli too)
                    % originally removed by because rotation wasn't working now kinda works with rotation, but sometimes doesn't
                    %% dosen't seem to create cl right?? get phase but no
                    %% sinusiod
                    imgSz = [imageRect(3)*2,imageRect(4)*2];
                    % Change clut for square gratings
                    if (get(handles.squaregratings,'Value')) % BA DOESTN WORK
                        [img cl] = generateSqGratings_lut(handles.orient(c),handles.freq(c),handles.TempFreq(c),handles.phase(c),handles.contrast(c),Duration, degPerPix,imgSz(1),imgSz(2),FrameHz,black,white,sizeLut);
                    else

                        [img cl] = generateGratings_lut(handles.orient(c),handles.freq(c),handles.TempFreq(c),handles.phase(c),handles.contrast(c),Duration, degPerPix,imgSz(1),imgSz(2),FrameHz,black,white,sizeLut);
                    end
                    destRect = windowRect;
                 else
                    [img cl] = generateCPGratings_lut(handles.orient(c),handles.freq(c),handles.TempFreq(c),handles.phase(c),handles.contrast(c),Duration, degPerPix,imageRect(3),imageRect(4),FrameHz,black,white,sizeLut);
                    condNAME = 'Counterphase Gratings';% BA

                end
                if exist('cl')
                    if c==1
                        cluts = zeros(256,3,size(cl,3),params.nCond);
                    end
                    cluts(:,:,:,c)=floor(cl);
                    textures(c,1)=Screen('MakeTexture',window,img);
                end
            end %cond
            fprintf('done generating')
            WAITTIME = 0; %BA
            FrameWait=1;

  
            %%%%% checkerboard  %%%%%%%%%%
        case 5
            condNAME = 'checkerboard';% BA

            params.nCond = size(handles.freq,2);
            textures = zeros(params.nCond,1);
            for c = 1:params.nCond
                [x y]= meshgrid(1:imageRect(3), 1:imageRect(4));
                contrast = handles.contrast(c);
                f= 2*pi*handles.freq(c)* degPerPix;
                img = 1+sign(sin(f*x).*sin(f*y));
                if contrast>1
                    contrast=1;
                end

                inc=(white-grey)*contrast;

                cl = offclut;
                cl(1,:) = grey-inc;
                cl(2,:) = grey+inc;
                cl(3,:) = grey+inc;

                fprintf('done generating')
                if c==1
                    cluts = zeros(256,3,size(cl,3),params.nCond);
                end
                cluts(:,:,:,c)=floor(cl);
                textures(c,1)=Screen('MakeTexture',window,img);
            end
            FrameWait = ceil(Duration*FrameHz);

            %%%%%%%  fullfield flash %%%%%%%%
        case 4 %
            condNAME = 'fullfield flash';% BA

            offclut(:,:)=black;
            textures = zeros(params.nCond,1);
            for c = 1:params.nCond
                img = ones(imageRect(4), imageRect(3));
                cl  = offclut;
                cl(:,:) = floor(white*handles.contrast(c));
                fprintf('done generating')
                if c==1
                    cluts = zeros(256,3,size(cl,3),params.nCond);
                end
                cluts(:,:,:,c)=floor(cl);
                textures(c,1)=Screen('MakeTexture',window,img);
            end % cond
            FrameWait = ceil(Duration*FrameHz);

            %%%%% drifting bars  %%%%%
        case 2
            condNAME = 'drifting bars';% BA

            textures = zeros(params.nCond,1);
            for c = 1:params.nCond

                % uncommented bva   040208
                frm = generateBars_blit(handles.orient(c),handles.freq(c),handles.speed(c),handles.contrast(c),handles.length(c), handles.positionX(c),Duration, degPerPix,imageRect(3),imageRect(4),FrameHz,black,white,2048);
                nFrames = size(frm,1);
                for f = 1:nFrames
                    textures(c,f)=Screen('MakeTexture',window,squeeze(frm(f,:,:)));
                end
                MovieRate = FrameHz;
                destRect = windowRect;
                clut=0;
                if c==1
                    save frames frm
                end
                %             end of uncomment bva 040208

                % %******* commented bva 040208
                %             [img cl] = generateBars_lut(handles.orient(c),handles.freq(c),handles.speed(c),handles.contrast(c),handles.length(c), handles.positionX(c),Duration, degPerPix,imageRect(3),imageRect(4),FrameHz,black,white,sizeLut);
                fprintf('done generating')
                if c==1
                    cluts = zeros(256,3,size(cl,3),params.nCond);
                end
                cluts(:,:,:,c)=floor(cl);
                textures(c,1)=Screen('MakeTexture',window,img);
                % %******* end of bva comment 040208

                % black background
                %offclut(:) = grey - (white-grey)*handles.contrast(c);
            end % cond
            FrameWait = 1;

            %%% flashing spots  %%%%%%%
        case 7
            condNAME = 'flashing spots';% BA

            params.nCond = size(handles.freq,2);
            textures = zeros(params.nCond,1);
            for c = 1:params.nCond
                [x y]= meshgrid(1:imageRect(3), 1:imageRect(4));
                contrast = handles.contrast(c);
                widthPix = handles.length(c)/degPerPix;
                posXpix = imageRect(3)/2 + handles.positionX(c)/degPerPix;
                posYpix = imageRect(4)/2 + handles.positionY(c)/degPerPix;

                img = double((x>(posXpix-widthPix/2)) & (x<(posXpix+widthPix/2)) & (y>(posYpix-widthPix/2)) & (y<(posYpix+widthPix/2)));
                if contrast>1
                    contrast=1;
                end
                inc = (white-grey)*contrast;
                cl = offclut;
                cl(2,:) = grey+inc;

                %%% non-grey background
                cl(1,:) = grey-0.75*inc;
                offclut(:) = grey - 0.75*inc;
                fprintf('done generating')

                if c==1
                    cluts = zeros(256,3,size(cl,3),params.nCond);
                end
                cluts(:,:,:,c)=floor(cl);
                textures(c,1)=Screen('MakeTexture',window,img);
            end % cond
            FrameWait = ceil(Duration*FrameHz);

            %%%%% movie %%%%%
        case 3
            load(get(handles.MovieName,'String'),'moviedata');
            condNAME = 'movie';% BA

            MovieMag = str2double(get(handles.MovieMag,'String'));
            MovieRate = str2double(get(handles.MovieRate,'String'));
            phasePeriod = str2double(get(handles.phasePeriod,'String'));

            % done loading movie

            length = str2double(get(handles.Length0,'String'));
            if length > 0
                length = length/MovieMag;
                length = length/degPerPix;
                moviedata = moviedata(1:round(length),:,:); %#ok
            end


            %         moviedata = moviedata(:,:,1:298); % BA just take a subset of the movie (%each part is 300 frames long)
            nFrames = size(moviedata,3);

            % fractional nReps truncates movie
            if (nReps < 1)
                nFrames = min(nFrames,floor(nFrames*nReps));
                nReps = 100000000;  %%%%% fix this
            end

            %% if multiple stimulus groups specified, break up movies into different stimulus conditions
            % period of the stimulus (i.e. one phase cycle)
            phasePeriodFrames = MovieRate * phasePeriod;
            if phasePeriodFrames == 0
                phasePeriodFrames = nFrames;
            end

            % number of different stimulus groups contained in movie, default 1
            % this does not constrain the number of repetitions in a group, which
            % is set by Duration and phasePeriodFrames
            stimgroups =str2double(get(handles.stimulusGroups,'String'));
            params.nCond = stimgroups;
            if params.nCond > 1
                % number of repetitions of each stimulus in one group
                nStimulusRepetitions = (nFrames/phasePeriodFrames) / params.nCond;
                disp(sprintf('%.1f conditions, %.1f stimuli per condition',params.nCond,nStimulusRepetitions));
                if (nStimulusRepetitions ~= floor(nStimulusRepetitions)) % must be integer
                    error('Movie not evenly divisible; is number of stimulus groups wrong?');
                end
            else
                params.nCond = size(handles.freq,2);   %%% if only one stimulus group, then use variables to set params.nCond
            end

            imageRect = SetRect(0,0,size(moviedata,1),size(moviedata,2));
            destRect = CenterRect(MovieMag*imageRect,windowRect);
            x0 = str2double(get(handles.PositionX0,'String'));
            if x0 ~= 0
                dx = x0/degPerPix;
                destRect = offsetrect(destRect,dx,0);
            end

            textures = zeros(1,nFrames);
            for f=1:nFrames
                textures(1,f)=Screen('MakeTexture',window,squeeze(moviedata(:,:,f))');
            end
            clear moviedata

            FrameWait = FrameHz/MovieRate;
            WAITTIME = 4; % BA this is a clug only relavent when clut=0 to make sure that DAQ side caputures the begining of the next stimulus presentation.
            %         % DAQ can miss trigger sometimes when acquiring very long, large
            %         files there is some time for saving etc.
    end


    % BA (setup keys and mouse
    KbName('UnifyKeyNames');
    UCkeys = declareUCkeys();

    % BA mask parameters TO DO ADD to GUI
    % mask radii
    % ADD to gui , and saved and reloaded
    UCparams.rx = str2num(get(handles.maskradiusx,'String')) ;% in pixels %% add change to degrs
    UCparams.ry = str2num(get(handles.maskradiusy,'String'));

    params.rStep = 5; % inc and descrement step mask radius, pixels
    params.RMAX = max(imageRect); % max radius of mask
    UCparams.lockMask = 1; % toggles mouse movement of mask

    %BA initial mask coordinates (later is mouse cursor coordinates)
    % ADD to gui and save and reloaded
    UCparams.mX =  str2num(get(handles.maskcenterx,'String')); % The x-coordinate of the mouse cursor
    UCparams.mY = str2num(get(handles.maskcentery,'String')); % The y-coordinate of the mouse cursor
    UCparams.rotation = 0;

    params.nMasks = 4; % number of cases in helperMakeMask
    UCparams.masktype = get(handles.popmenuMask,'Value')-1;
    UCparams.masktex =  helperMakeMask(UCparams.rx,UCparams.ry,window,UCparams.masktype);


    % ADD to GUI control of this
    UCparams.bautoChangeContrast = 1;

    params.window = window; % for passing into keycontrol function


    %%%% clear screen
    Screen('FillRect',window,128);
    Screen('DrawText',window,sprintf('Finished stimuli'),10,40);
    Screen('Flip',window);

    %%%% gamma correction
    flat_clut = [(0:1/255:1)' (0:1/255:1)' (0:1/255:1)'];
    if ~isempty (handles.USER.GAMMATABLE); load(handles.USER.GAMMATABLE )
        screen('LoadNormalizedGammaTable',window,inv_gamma_clut);
        clear spyderCaldata
    else % correct with out calibration
        screen_gamma=2; % commented by BA on 062208 (although this calibration is
        % % pretty` close to right(within a few percent)
        gamma_clut = flat_clut.^(1/screen_gamma);
%         screen('LoadNormalizedGammaTable',window,gamma_clut);
        screen('LoadNormalizedGammaTable',window,flat_clut);
    end
    %%% get number of frames
    if clut
        nFrames = size(cluts,3);
    else
        nFrames = size(textures,2);
    end

    %% setup synchronization

    statusfile = fopen('statusfile.txt','w');
    startTime = GetSecs();
 
     sync = handles.USER.Sync.ktdtSync; % sync is a left over from Cris's code where there were different sync sources. BA

    if sync == handles.USER.Sync.ktdtSync
        dio = digitalio('parallel','LPT1');
        hwline = addline(dio,0:1,2,'out','bitLine');    %% pins 1,14
        addline(dio,0:7,'out'); %% pins 2-9

        %%% for some reason, the slowest part of putvalue when using parallel port
        %%% is this step, finding the parent uddobj.
        %%% So we look this up in the beginning, and then directly call the
        %%% putvalue function with uddobject, data, and line numbers.
        %%% Not exactly sure why this works, but it reduces time per call
        %%% from 2msec to 20usec (at least in previous versions of daqtoolbox
        parent = get(dio.bitLine, 'Parent');
        parentuddobj = daqgetfield(parent{1},'uddobject');
        % stimLine = 1;
        % frameLine = 2;
        bitLine=1:2;
        condNum=3:10;
        bitOn=0;
        bitOff=1;
        putvalue(parentuddobj,0,condNum);
        putvalue(parentuddobj,[bitOff bitOff],bitLine); %
    end




    %% set background clut
    if clut
        moglClutBlit(window,textures(1),offclut);
        % currentclut=offclut;
    end

    clearBkgrnd = get(handles.bkgrnd,'Value');

    %% add blank stimulus as extra condition
    if get(handles.blankstim,'Value')
        if clut
            params.nCond = params.nCond + 1;
            cluts(:,:,:,params.nCond) = offclut(1,1);
            textures(params.nCond) = textures(params.nCond-1);
        else  %%% no need for blank in movies (at least for now)
            sprintf('no blank available for movies')
        end
    end

    %% add full field flicker as extra condition, for drifiting gratings
    if get(handles.FullFlicker,'Value') && (get(handles.StimType,'Value')==1)
        params.nCond = params.nCond + 1;
        if clut
            cluts(:,:,:,params.nCond)=cluts(:,:,:,1);
            textures(params.nCond)=Screen('MakeTexture',window,ones(size(img)));
        end
    end

    %% set up run variables

    FrameInt = 1/FrameHz;
    WaitInt = str2double(get(handles.WaitInterval,'String'));

    s1 = zeros(nFrames,1);
    ds = zeros(nFrames-1,1);

    if stim == 3  % movie condition
        vblstore = nan(1,min(nFrames*params.nCond*20,1000000)); % BA predefine space for saving vbl time (this is a backup in case
        vbli = 1;
    else
        vblstore = -1;
    end

    %% BA for saving condList
    MAXSHUFFLES = 500; %number of times that condList will be reshuffled, after which the first shuffle will be used again
    condListStore = nan(1,params.nCond*MAXSHUFFLES);
    % texturec = 1;

    %% finally, run stimulus!!
    warning off MATLAB:concatenation:integerInteraction  %%%% this error comes up in generating udp packet
    %     ListenChar(2); % BA removed this cause it is annoying

    % Make sure all Rushed variables and functions are in memory
    % before raising priority
    UCparams.doneStim = 0;
    UCparams.break = 0;
    iter = 0;
    numiters = nReps * params.nCond * nStimulusRepetitions;
    % stimulusrep applies when there are multiple stimuli for a particular
    % condition, e.g. different noise patterns.
    stimulusrep = 1;

    GetSecs;
    Screen('Screens');
    %     HideCursor; %BA

    lastsecs = [];



    %%% loop on conditions
    while ~UCparams.doneStim
        %%% randomize conditions
        if get(handles.randomize,'Value');
            if mod(iter,params.nCond)==0   %%% shuffle condition list each repeat
                if iter>0 && mod(iter,params.nCond*MAXSHUFFLES)==0 %number of times that condList has been reshuffled exceeds MAXSHUFFLES resuse previous REShuffles in order
                    temp =  mod(iter,params.nCond*MAXSHUFFLES);
                    temp = temp + ~temp; % make start from 1
                    condList = condListStore(temp:temp+params.nCond+1); %
                else
                    condList = Shuffle(1:params.nCond);
                    condListStore(max(iter,1):max(iter,1)+params.nCond-1) = condList; % save so that is saved to disk
                end
            end
        else
            condList = 1:params.nCond;
            condListStore(1:params.nCond) = condList;
        end

        %%% choose condition for this iteration and send it out
        if UCparams.bautoChangeContrast % BA allows manual control of condition
            UCparams.c = condList(mod(iter,params.nCond)+1);
        end
        if sync == handles.USER.Sync.ktdtSync
            putvalue(parentuddobj, UCparams.c,condNum);
        end
        %

        %

        %%% raise the priority
        priorityLevel = MaxPriority(window);
        Priority(priorityLevel);
        %%% loop for blit anim              ation
        if ~clut
            if stim==3 % BA not sure what this is for don't like it for my gratings but leaving it for moviedata
                vbl = Screen('Flip',window);  %%% initial flip, to sync with vertical blank
            end
            %% set which frames to show
            if (stim == 3 && stimgroups > 1)
                % if multiple condition movie, randomly choose a condition
                offset = ( UCparams.c-1)*nStimulusRepetitions*phasePeriodFrames + ...
                    (stimulusrep-1)*phasePeriodFrames;
                minframe = offset + 1;
                maxframe = offset + phasePeriodFrames;
            else
                % show movie all the way through regardless of randomization
                minframe = 1;
                maxframe = nFrames;
            end


            fwrite(u,[fname '*' condNAME]);% BA output stimulus data to DAQ PC
            %             putvalue(parentuddobj,[bitOn bitOff],bitLine); % BA add delay after movie starts so that it is obvious looking at the frames data that the movie just started
            WaitSecs(WAITTIME);% BA (delay so that next trigger is caught by DAQ) kluggy, WAITTIME>0 only for reverse correlation


            %%% loop through framess
            for f = minframe:maxframe
                if UCparams.break
                    break;
                end
                s1(f) = GetSecs;
                Screen('DrawTexture',window, textures(UCparams.c,f),[],destRect,UCparams.rotation);
                Screen('DrawTexture', window, UCparams.masktex, [],CenterRectOnPoint(windowRect*2, UCparams.mX, UCparams.mY));                         % draw mask
                %% BA CHECK may be drawing mask outside of window
                Screen('DrawText',window,sprintf('C%d %d %d',params.nCond,iter,stimulusrep),10,2,grey);
                frameChangeBox(window,windowRect); % adds box to corner of screen and switch box from gray to white everytime function is called 

                putvalue(parentuddobj,[bitOn bitOn],bitLine); 
                if f > 1
                    vbl = Screen('Flip',window, vbl + (FrameWait - 0.5) * FrameInt);
                    % BA why minus 0.5 (seems arbitrary)
                else
                    vbl = Screen('Flip',window);
                end
                putvalue(parentuddobj,[bitOn bitOff],bitLine);
                vblstore(vbli)= vbl;

                % BA USer Control stimulus and mask with keyboard and mouse
                if ~UCparams.lockMask
                    [UCparams.mX, UCparams.mY, UCparams.buttons] = GetMouse;
                end

                [keyIsDown, secs, keyCode] = KbCheck;
                if isempty(lastsecs)|( secs -lastsecs) >= 1/8; % set a maximum rate at which a key can be pressed otherwise holding button makes bar spin to fast
                    if keyIsDown %%% charavail would be much better, but doesn't seem to work
                        lastsecs= secs;
                        UCparams = helperUserControl(keyCode,UCkeys,UCparams,params);
                    end
                end


            end

            %%% done with stimulus
            if clearBkgrnd
                Screen('FillRect',window,grey);
            else % BA
                if stim==3 % BA not sure what this is for don't like it for my gratings but leaving it for moviedata

                    Screen('DrawTexture',window, textures(UCparams.c,f),[],destRect,UCparams.rotation);
                    Screen('DrawTexture', window, UCparams.masktex, [],CenterRectOnPoint(windowRect*2, UCparams.mX, UCparams.mY));
                end
            end

            if stim==3 % BA not sure what this is for don't like it for my gratings but leaving it for moviedata
                frameChangeBox(window,windowRect); % adds box to corner of screen and switch box from gray to white everytime function is called
                vbl = Screen('Flip',window, vbl+ (FrameWait - 0.5) * FrameInt);
            end
            if sync == handles.USER.Sync.ktdtSync
                putvalue(parentuddobj,[bitOff bitOff],bitLine);
            end

            %%% loop for lookup table animation
        elseif clut
            clutcond = (squeeze(cluts(:,:,:, UCparams.c)));
            %%% first clut loaded is slow, so must load something
            %%% (at least in old version)
            %             moglClutBlit(window,textures( UCparams.c,1),currentclut);
            %             vbl= Screen('Flip',window);
            fwrite(u,[fname '*' condNAME]);% BA output stimulus data to DAQ PC

            %%%% loop through frames
            for f = 1:nFrames
                if UCparams.break
                    break;
                end
                s1(f) = GetSecs;
                %                moglClutBlitBA(window,textures( UCparams.c),clutcond(:,:,f),UCparams.rotation);
                moglClutBlit(window,textures( UCparams.c),clutcond(:,:,f),UCparams.rotation);
                %                moglClutBlit(window,textures(UCparams.c),clutcond(:,:,f));
                Screen('DrawText',window,sprintf('%dd %1.2fcpd %1.1fhz %1.2f C%d %d %d',handles.orient(UCparams.c)+UCparams.rotation,handles.freq(UCparams.c),handles.TempFreq(UCparams.c),handles.contrast(UCparams.c),params.nCond,iter,stimulusrep),10,2,grey);

                % draw mask
                Screen('DrawTexture', window, UCparams.masktex, [],CenterRectOnPoint(windowRect*2, UCparams.mX, UCparams.mY));
                frameChangeBox(window,windowRect); % adds box to corner of screen and switch box from gray to white everytime function is called 

                putvalue(parentuddobj,[bitOn bitOn],bitLine); % BA technically sync == handles.USER.Sync.ktdtSync only
                if f > 1
                    % BA TODO send parrellel port bits around this for
                    % better measure of frame
                    vbl = Screen('Flip',window, vbl + (FrameWait - 0.5) * FrameInt);
                else
                    vbl = Screen('Flip',window);
                end
                putvalue(parentuddobj,[bitOn bitOff],bitLine);

                % BA USer Control stimulus and mask with keyboard and mouse
                if ~UCparams.lockMask
                    [UCparams.mX, UCparams.mY, UCparams.buttons] = GetMouse;
                end

                [keyIsDown, secs, keyCode] = KbCheck;
                if isempty(lastsecs)|( secs -lastsecs) >= 1/10; % set a maximum rate at which a key can be pressed otherwise holding button makes bar spin to fast
                    if keyIsDown %%% charavail would be much better, but doesn't seem to work
                        lastsecs= secs;
                        UCparams = helperUserControl(keyCode,UCkeys,UCparams,params,handles);
                    end
                end

            end
            %%% done with stimulus

            if clearBkgrnd
                moglClutBlit(window,textures( UCparams.c),offclut);
                % currentclut = offclut;
            else % BA why not do NOTHING instead?
                %                 moglClutBlitBA(window,textures( UCparams.c),clutcond(:,:,nFrames),UCparams.rotation);
                %                 moglClutBlit(window,textures( UCparams.c),clutcond(:,:,nFrames));
                moglClutBlit(window,textures( UCparams.c),clutcond(:,:,nFrames),UCparams.rotation);
                Screen('DrawTexture', window, UCparams.masktex, [],CenterRectOnPoint(windowRect*2, UCparams.mX, UCparams.mY)); % mask BA
                Screen('DrawTexture', window, UCparams.masktex, [],CenterRectOnPoint(windowRect*2, UCparams.mX, UCparams.mY)); % mask BA
                % currentclut = squeeze(clutcond(:,:,nFrames));
            end
            frameChangeBox(window,windowRect); % adds box to corner of screen and switch box from gray to white everytime function is called

            vbl = Screen('Flip',window, vbl + (FrameWait - 0.5) * FrameInt);

            if sync == handles.USER.Sync.ktdtSync
                putvalue(parentuddobj,[bitOff bitOff],bitLine);
            end
        end    %%%clut

        Priority(0);
        WaitSecs(WaitInt);

        if nFrames > 1
            ds = max(ds,diff(s1));
        end

        iter = iter + 1;

        elapsedTime = GetSecs - startTime;
        fprintf(statusfile,'%d %d %0.2f \r\n',iter,int16(numiters),elapsedTime);

        if mod(iter,params.nCond) == 0 %% done all conditions, move on to next stimulus
            stimulusrep = stimulusrep + 1;
            % if new Rep of entire movie, start new stimulusrep
            stimulusrep = mod(stimulusrep-1,nStimulusRepetitions)+1;
        end

        % look for stop message on keyboard
        % BA commented out New key entries are handled in frame loop
        %         [keyIsDown, secs, keyCode] = KbCheck;
        %         if keyIsDown %%% charavail would be much better, but doesn't seem to work
        %             UCparams.doneStim = 1;
        %             keyspressed = KbName(find(keyCode));
        %             disp(sprintf('Exit on %s key pressed',keyspressed{1}));
        %         end


        % test for all stimuli complete
        if (iter >= numiters)
            UCparams.doneStim = 1;
            disp('Exit on completion');
        end

    end %while ~doneStim
    Priority(0);

    %%% Save conditions and condition order to disk
    condListStore = condListStore(~isnan(condListStore));
    save([fname '_Condstore'],'condListStore','iter'); % this is a back up (each condition should be transmitted to DAQ PC each iteration)
    %%% Save vbl data to disk %%%

    vblstore = vblstore(~isnan(vblstore));
    save([fname '_vblstore'],'vblstore');

    %%% UPDATE GUI values BA%%%%%%%%%%
    set(handles.popmenuMask,'Value',UCparams.masktype+1);
    set(handles.maskradiusx,'String',num2str(UCparams.rx)); % in pixels
    set(handles.maskradiusy,'String',num2str(UCparams.ry));
    set(handles.maskcenterx,'String',num2str(UCparams.mX)); % The x-coordinate of the mouse cursor
    set(handles.maskcentery,'String',num2str(UCparams.mY)); % The y-coordinate of the mouse cursor



    %%%% cleanup   %%%%
    if exist('u','var');        fclose(u); delete(u); clear u  ;  end
    moglClutBlit;
    ListenChar(1);
    pnet('closeall')
    fclose(statusfile);
    Screen('LoadNormalizedGammaTable',window,flat_clut);
    screen('CloseAll');

    f = figure(1);
    set(f,'Position',[1 35 1024 130])
    plot(ds);
    title('Dropped frames');

    %BA save frames timing
    ShowCursor;

    %%% if there's an error, clean up and rethrow the error
catch ME
    getReport(ME)
    
    if exist('u','var');        fclose(u); delete(u);   clear u ;   end
    moglClutBlit;
    ShowCursor;
    ListenChar(1);

    pnet('closeall')
    if exist('statusfile','var')
        fclose(statusfile)
    end
    if exist('flat_clut','var')
        Screen('LoadNormalizedGammaTable',window,flat_clut);
    end
    if clut
        moglClutBlit;  %%%% need to close this, or won't work next time
    end
    Screen('CloseAll');

    Priority(0);
    psychrethrow(psychlasterror);
end

function frameChangeBox(window,imageRect)
% function to add box to corner of screen. box switches everytime fucntion
% is called. use to monitor frame changes. ie. call just before each flip, thus provides an location where a photodiode can monitor the
% frame changes (EF's idea)
persistent c;
if isempty(c)||c==255;
    c = 128;
else
    c = 255;
end
Screen('FillRect',window,c,[imageRect(3)-50 imageRect(4)-50 imageRect(3) imageRect(4)]);

% --- Executes on button press in SaveParams.
function SaveParams_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to SaveParams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fname, pname] = uiputfile('*.mat','Parameter File',handles.USER.paramdirpath);
if (fname == 0), return; end % canceled

fname = fullfile(pname,fname);
SaveParams(handles,fname);

handles.USER.paramdirpath = pname;

%--- function to save parameters, called by SaveParams, or on Run_Btn
function SaveParams(handles,fname)

Orient0 = str2double(get(handles.Orient0,'String')); %#ok
Freq0 = str2double(get(handles.Freq0,'String')); %#ok
Speed0 = str2double(get(handles.Speed0,'String')); %#ok
Contrast0 = str2double(get(handles.Contrast0,'String')); %#ok
TempFreq0 = str2double(get(handles.TempFreq0,'String')); %#ok
Duration= str2double(get(handles.Duration,'String')); %#ok
Phase0 = str2double(get(handles.Phase0,'String')); %#ok
Length0 = str2double(get(handles.Length0,'String')); %#ok
PositionX0 = str2double(get(handles.PositionX0,'String')); %#ok
PositionY0 = str2double(get(handles.PositionY0,'String')); %#ok
WaitInterval = str2double(get(handles.WaitInterval,'String')); %#ok
if isfield(handles,'eyeCond0'); eyeCond0 = str2double(get(handles.eyeCond0,'String')); else eyeCond0 = 0; end; % removed eye (left this for backward compatibility)

StimulusStr = get(handles.StimType,'String'); %#ok
StimulusNum = get(handles.StimType,'Value'); %#ok

PixelsX = str2double(get(handles.PixelsX,'String')); %#ok
PixelsY = str2double(get(handles.PixelsY,'String')); %#ok
SizeX = str2double(get(handles.SizeX,'String')); %#ok
SizeY = str2double(get(handles.SizeY,'String')); %#ok
ScreenDist = str2double(get(handles.ScreenDist,'String')); %#ok

Var1Str = get(handles.Var1,'String'); %#ok
Var1Val = get(handles.Var1,'Value'); %#ok

Start1 = str2double(get(handles.Start1,'String')); %#ok
Stop1 = str2double(get(handles.Stop1,'String')); %#ok
nSteps1 = str2double(get(handles.nSteps1,'String')); %#ok
LinLog1 = get(handles.LinLog1,'Value'); %#ok


Var2Str = get(handles.Var2,'String'); %#ok
Var2Val = get(handles.Var2,'Value'); %#ok

Start2 = str2double(get(handles.Start2,'String')); %#ok
Stop2 = str2double(get(handles.Stop2,'String')); %#ok
nSteps2 = str2double(get(handles.nSteps2,'String')); %#ok
LinLog2 = get(handles.LinLog2,'Value'); %#ok

MovieName = get(handles.MovieName,'String'); %#ok
MovieMag = str2double(get(handles.MovieMag,'String')); %#ok
MovieRate = str2double(get(handles.MovieRate,'String')); %#ok
phasePeriod = str2double(get(handles.phasePeriod,'String')); %#ok
stimulusGroups = str2double(get(handles.stimulusGroups,'String')); %#ok

if isfield(handles,'SyncSource');SyncSource = get(handles.SyncSource,'Value');else SyncSource = 0; end; % removed SyncSource (left this for backward compatibility)
 %#ok

blankbkgrnd = get(handles.bkgrnd,'Value'); %#ok
randomize = get(handles.randomize,'Value'); %#ok
blankstim = get(handles.blankstim,'Value'); %#ok
FullFlicker = get(handles.FullFlicker,'Value'); %#ok
nReps = get(handles.nReps,'String'); %#ok

orient = handles.orient; %#ok
spfreq = handles.freq; %#ok
speed = handles.speed; %#ok
contrast = handles.contrast; %#ok
phase = handles.phase; %#ok
TempFreq = handles.TempFreq; %#ok
positionX = handles.positionX; %#ok
positionY = handles.positionY; %#ok
length = handles.length; %#ok
squaregratings = get(handles.squaregratings,'Value'); %#ok

% mask BA
maskstr = get(handles.popmenuMask,'String'); %#ok
popmenuMask = get(handles.popmenuMask,'Value');
maskcenterx =  get(handles.maskcenterx,'String');
maskcentery = get(handles.maskcentery,'String');
maskcenterdeg = get(handles.maskcenterdeg,'String');
maskradiusx = get(handles.maskradiusx,'String');
maskradiusy =get(handles.maskradiusy,'String');
maskmeanradiusdeg = get(handles.maskmeanradiusdeg,'String');


save(fname, 'Orient0', 'Freq0', 'TempFreq0','Phase0','Speed0', 'Contrast0','Duration','Length0','PositionX0', 'PositionY0', 'WaitInterval','eyeCond0','StimulusStr', 'StimulusNum', ...
    'PixelsX','PixelsY','SizeX','SizeY','ScreenDist','Var1Str','Var1Val','Start1','Stop1','nSteps1','LinLog1', ...
    'Var2Str','Var2Val','Start2','Stop2','nSteps2','LinLog2','MovieName','MovieMag','MovieRate', ...
    'orient', 'spfreq', 'speed', 'contrast','phase','TempFreq','length','positionX','positionY','blankbkgrnd','randomize','blankstim','FullFlicker',...
    'nReps','phasePeriod','stimulusGroups','squaregratings','maskstr','popmenuMask','maskcenterx','maskcentery','maskcenterdeg','maskradiusx','maskradiusy','maskmeanradiusdeg');

% --- Executes on button press in LoadParams.
function LoadParams_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to LoadParams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fname, pname] = uigetfile('*.mat','Parameter File',handles.USER.paramdirpath);
if (fname == 0), return; end % canceled

load(fullfile(pname,fname));
handles.USER.paramdirpath = pname;

set(handles.StimType,'Value',StimulusNum);
set(handles.Orient0,'String',num2str(Orient0));
set(handles.Freq0,'String',num2str(Freq0));
set(handles.Speed0,'String',num2str(Speed0));
set(handles.Contrast0,'String',num2str(Contrast0));
set(handles.Duration,'String',num2str(Duration));

set(handles.PixelsX,'String',num2str(PixelsX));
set(handles.PixelsY,'String',num2str(PixelsY));
set(handles.SizeX,'String',num2str(SizeX));
set(handles.SizeY,'String',num2str(SizeY));
set(handles.ScreenDist,'String',num2str(ScreenDist));

set(handles.Var1,'Value',Var1Val);
set(handles.Start1,'String',num2str(Start1));
set(handles.Stop1,'String',num2str(Stop1));
set(handles.nSteps1,'String',num2str(nSteps1));
set(handles.LinLog1,'Value',LinLog1);

set(handles.Var2,'Value',Var2Val);
set(handles.Start2,'String',num2str(Start2));
set(handles.Stop2,'String',num2str(Stop2));
set(handles.nSteps2,'String',num2str(nSteps2));
set(handles.LinLog2,'Value',LinLog2);

set(handles.MovieName,'String',MovieName);
set(handles.MovieMag,'String',num2str(MovieMag));
set(handles.MovieRate,'String',num2str(MovieRate));

if exist('TempFreq0','var')   %added 28Mar2006
    set(handles.TempFreq0,'String',num2str(TempFreq0));
    set(handles.Phase0,'String',num2str(Phase0));
end

if exist('Length0','var')  %% added 20Jul2006
    set(handles.randomize,'Value',randomize);
    set(handles.blankstim,'Value',blankstim);
    set(handles.bkgrnd,'Value',blankbkgrnd);
    set(handles.Length0,'String',num2str(Length0));
    %  set(handles.PositionX0,'String',num2str(Position0));
    set(handles.WaitInterval,'String',num2str(WaitInterval));
end

if exist('PositionY0','var') %%% added 18Aug2006
    set(handles.PositionX0,'String',num2str(PositionX0));
    set(handles.PositionY0,'String',num2str(PositionY0));
end

if exist('FullFlicker','var') %%% added Jan172007
    set(handles.FullFlicker,'Value',FullFlicker);
end

if exist('nReps','var') %%% added 15Aug2007
    set(handles.nReps,'String',num2str(nReps));
    set(handles.phasePeriod,'String',num2str(phasePeriod));
    set(handles.stimulusGroups,'String',num2str(stimulusGroups));
end

if exist('squaregratings','var') %%% added 27Sep2007
    set(handles.squaregratings,'Value',squaregratings);
end

% if exist('eyeCond0','var') %%% added 10Oct2007
%     set(handles.eyeCond0,'String',num2str(eyeCond0));
% end

if exist('SyncSource','var') %%% added 20Oct2007
    set(handles.squaregratings,'Value',SyncSource);
end

if exist('popmenuMask','var') % %BA 081909
    set(handles.popmenuMask,'Value',popmenuMask)
    set(handles.maskcenterx,'String',maskcenterx)
    set(handles.maskcentery,'String',maskcentery)
    set(handles.maskcenterdeg,'String',maskcenterdeg)
    set(handles.maskradiusx,'String',maskradiusx)
    set(handles.maskradiusy,'String',maskradiusy)
    set(handles.maskmeanradiusdeg,'String',maskmeanradiusdeg)
end


StimType_Callback(handles.StimType,eventdata,handles);
Var1_Callback(handles.Var1,eventdata,handles);
Var2_Callback(handles.Var2,eventdata,handles);


% --- Executes on button press in WaitSync.
function WaitSync_Callback(hObject, eventdata, handles) %#ok



function TempFreq0_Callback(hObject, eventdata, handles) %#ok
StimType_Callback(handles.StimType,eventdata,handles);

% --- Executes during object creation, after setting all properties.
function TempFreq0_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Phase0_Callback(hObject, eventdata, handles) %#ok
StimType_Callback(handles.StimType,eventdata,handles);

% --- Executes during object creation, after setting all properties.
function Phase0_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Orient0_Callback(hObject, eventdata, handles) %#ok
StimType_Callback(handles.StimType,eventdata,handles);

% --- Executes during object creation, after setting all properties.
function Orient0_CreateFcn(hObject, eventdata, handles) %#ok
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function Freq0_Callback(hObject, eventdata, handles) %#ok
StimType_Callback(handles.StimType,eventdata,handles);

% --- Executes during object creation, after setting all properties.
function Freq0_CreateFcn(hObject, eventdata, handles) %#ok
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function Speed0_Callback(hObject, eventdata, handles) %#ok
StimType_Callback(handles.StimType,eventdata,handles);

% --- Executes during object creation, after setting all properties.
function Speed0_CreateFcn(hObject, eventdata, handles) %#ok
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function Contrast0_Callback(hObject, eventdata, handles) %#ok
StimType_Callback(handles.StimType,eventdata,handles);

% --- Executes during object creation, after setting all properties.
function Contrast0_CreateFcn(hObject, eventdata, handles) %#ok
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function Duration_Callback(hObject, eventdata, handles) %#ok

function Duration_CreateFcn(hObject, eventdata, handles) %#ok
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



% --- Executes on selection change in Var1.
function Var1_Callback(hObject, eventdata, handles)
% hObject    handle to Var1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Var1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Var1

if get(hObject,'Value')==1 %none
    set(handles.Start1,'Enable','off');
    set(handles.Stop1,'Enable','off');
    set(handles.nSteps1,'Enable','off');
    set(handles.LinLog1,'Enable','off');
    set(handles.Var1Range,'Enable','off');
else
    set(handles.Start1,'Enable','on');
    set(handles.Stop1,'Enable','on');
    set(handles.nSteps1,'Enable','on');
    set(handles.LinLog1,'Enable','on');
    set(handles.Var1Range,'Enable','on');
end

Var1Range_Callback(handles.Var1Range,eventdata,handles);


% --- Executes during object creation, after setting all properties.
function Var1_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to Var1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in Var2.
function Var2_Callback(hObject, eventdata, handles)
% hObject    handle to Var2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Var2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Var2

if get(hObject,'Value')==1 %none
    set(handles.Start2,'Enable','off');
    set(handles.Stop2,'Enable','off');
    set(handles.nSteps2,'Enable','off');
    set(handles.LinLog2,'Enable','off');
    set(handles.Var2Range,'Enable','off');

else
    set(handles.Start2,'Enable','on');
    set(handles.Stop2,'Enable','on');
    set(handles.nSteps2,'Enable','on');
    set(handles.LinLog2,'Enable','on');
    set(handles.Var2Range,'Enable','on');
end

Var2Range_Callback(handles.Var2Range,eventdata,handles);

% --- Executes during object creation, after setting all properties.
function Var2_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to Var2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in LinLog1.
function LinLog1_Callback(hObject, eventdata, handles) %#ok
Var1_Callback(handles.Var1, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function LinLog1_CreateFcn(hObject, eventdata, handles) %#ok
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function Start1_Callback(hObject, eventdata, handles) %#ok

Var1_Callback(handles.Var1, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function Start1_CreateFcn(hObject, eventdata, handles) %#ok
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function Stop1_Callback(hObject, eventdata, handles) %#ok
Var1_Callback(handles.Var1, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function Stop1_CreateFcn(hObject, eventdata, handles) %#ok
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function nSteps1_Callback(hObject, eventdata, handles) %#ok
Var1_Callback(handles.Var1, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function nSteps1_CreateFcn(hObject, eventdata, handles) %#ok
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function Start2_Callback(hObject, eventdata, handles) %#ok
Var2_Callback(handles.Var2, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function Start2_CreateFcn(hObject, eventdata, handles) %#ok
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function Stop2_Callback(hObject, eventdata, handles) %#ok
Var2_Callback(handles.Var2, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function Stop2_CreateFcn(hObject, eventdata, handles) %#ok
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function nSteps2_Callback(hObject, eventdata, handles) %#ok
Var2_Callback(handles.Var2, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function nSteps2_CreateFcn(hObject, eventdata, handles) %#ok
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in LinLog2.
function LinLog2_Callback(hObject, eventdata, handles) %#ok
Var2_Callback(handles.Var2, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function LinLog2_CreateFcn(hObject, eventdata, handles) %#ok
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function Var1Range_Callback(hObject, eventdata, handles) %#ok

nSteps1 = str2double(get(handles.nSteps1,'String'));
Start1 = str2double(get(handles.Start1,'String'));
Stop1 = str2double(get(handles.Stop1,'String'));

variscircular = 0;
if (get(handles.Var1,'Value')) == 2 && (Start1 == 0) && (Stop1 == 360)
    variscircular = 1; % orientation
end

if get(handles.LinLog1,'Value')==1
    if variscircular
        Var1Range = linspace(Start1, Stop1, nSteps1+1);
        Var1Range = Var1Range(1:end-1);
    else
        Var1Range = linspace(Start1, Stop1, nSteps1);
    end
else
    Var1Range = logspace(log10(Start1), log10(Stop1), nSteps1);
end
set(hObject,'String',mat2str(Var1Range,3));

[handles.orient handles.freq handles.speed handles.contrast handles.phase ...
    handles.TempFreq handles.var1value handles.var2value handles.positionX handles.positionY handles.length handles.eye] = generateVarParams(handles);

codes = [0 4 23 10 13 24 18 7 8 22];
handles.var1code = codes(get(handles.Var1,'Value'));
handles.var2code = codes(get(handles.Var2,'Value'));

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function Var1Range_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to Var1Range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function Var2Range_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to Var2Range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Var2Range as text
%        str2double(get(hObject,'String')) returns contents of Var2Range as a double

nSteps2 = str2double(get(handles.nSteps2,'String'));
Start2 = str2double(get(handles.Start2,'String'));
Stop2 = str2double(get(handles.Stop2,'String'));

variscircular = 0;
if (get(handles.Var2,'Value')) == 2 && (Start2 == 0) && (Stop2 == 360)
    variscircular = 1; % orientation
end

if get(handles.LinLog2,'Value')==1
    if variscircular
        Var2Range = linspace(Start2, Stop2, nSteps2+1);
        Var2Range = Var2Range(1:end-1);
    else
        Var2Range = linspace(Start2, Stop2, nSteps2);
    end
else
    Var2Range = logspace(log10(Start2), log10(Stop2), nSteps2);
end
set(hObject,'String',mat2str(Var2Range,3));

[handles.orient handles.freq handles.speed handles.contrast handles.phase ...
    handles.TempFreq handles.var1value handles.var2value handles.positionX handles.positionY handles.length handles.eye] = generateVarParams(handles);

codes = [0 4 23 10 13 24 18 7 8  22];
handles.var1code =codes(get(handles.Var1,'Value'));
handles.var2code =codes(get(handles.Var2,'Value'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function Var2Range_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to Var2Range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function MovieName_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to MovieName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MovieName as text
%        str2double(get(hObject,'String')) returns contents of MovieName as a double

updateparamsfrommovie(get(hObject,'String'),handles);

% --- Executes during object creation, after setting all properties.
function MovieName_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to MovieName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in SelectMovieName.
function SelectMovieName_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to SelectMovieName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


[fname pname] = uigetfile('.mat','Movie Name',handles.USER.moviedirpath);
if (fname == 0); return; end

fullmoviename = fullfile(pname,fname);
handles.USER.moviedirpath = pname;

set(handles.MovieName,'String',fullmoviename);

updateparamsfrommovie(fullmoviename,handles);


function MovieMag_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to MovieMag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MovieMag as text
%        str2double(get(hObject,'String')) returns contents of MovieMag as a double


% --- Executes during object creation, after setting all properties.
function MovieMag_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to MovieMag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in SyncSource.
function SyncSource_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to SyncSource (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns SyncSource contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SyncSource

% --- Executes during object creation, after setting all properties.
function SyncSource_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to SyncSource (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function WaitInterval_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to WaitInterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WaitInterval as text
%        str2double(get(hObject,'String')) returns contents of WaitInterval as a double


% --- Executes during object creation, after setting all properties.
function WaitInterval_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to WaitInterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function ScreenNum_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to ScreenNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ScreenNum as text
%        str2double(get(hObject,'String')) returns contents of ScreenNum as a double
FrameHz_Callback(handles.FrameHz,eventdata,handles);
PixelsX_Callback(handles.PixelsX,eventdata,handles);
PixelsY_Callback(handles.PixelsY,eventdata,handles);



% --- Executes during object creation, after setting all properties.
function ScreenNum_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to ScreenNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




function PixelsX_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to PixelsX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PixelsX as text
%        str2double(get(hObject,'String')) returns contents of PixelsX as a double
screenNum = str2double(get(handles.ScreenNum,'String'));
rect = Screen(screenNum,'Rect');
set(hObject,'String',num2str(rect(3)));


% --- Executes during object creation, after setting all properties.
function PixelsX_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to PixelsX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function PixelsY_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to PixelsY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PixelsY as text
%        str2double(get(hObject,'String')) returns contents of PixelsY as a double
screenNum = str2double(get(handles.ScreenNum,'String'));
rect = Screen(screenNum,'Rect');
set(hObject,'String',num2str(rect(4)));

% --- Executes during object creation, after setting all properties.
function PixelsY_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to PixelsY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function SizeX_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to SizeX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SizeX as text
%        str2double(get(hObject,'String')) returns contents of SizeX as a double


% --- Executes during object creation, after setting all properties.
function SizeX_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to SizeX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function SizeY_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to SizeY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SizeY as text
%        str2double(get(hObject,'String')) returns contents of SizeY as a double


% --- Executes during object creation, after setting all properties.
function SizeY_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to SizeY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function ScreenDist_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to ScreenDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ScreenDist as text
%        str2double(get(hObject,'String')) returns contents of ScreenDist as a double


% --- Executes during object creation, after setting all properties.
function ScreenDist_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to ScreenDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function FrameHz_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to FrameHz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FrameHz as text
%        str2double(get(hObject,'String')) returns contents of FrameHz as a double
set(hObject,'String',num2str(FrameRate(str2double(get(handles.ScreenNum,'String')))));


% --- Executes during object creation, after setting all properties.
function FrameHz_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to FrameHz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function MovieRate_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to MovieRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MovieRate as text
%        str2double(get(hObject,'String')) returns contents of MovieRate as a double


% --- Executes during object creation, after setting all properties.
function MovieRate_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to MovieRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in bkgrnd.
function bkgrnd_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to bkgrnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bkgrnd


% --- Executes on button press in randomize.
function randomize_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to randomize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of randomize


function Length0_Callback(hObject, eventdata, handles) %#ok
StimType_Callback(handles.StimType,eventdata,handles);

% --- Executes during object creation, after setting all properties.
function Length0_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function PositionX0_Callback(hObject, eventdata, handles) %#ok
StimType_Callback(handles.StimType,eventdata,handles);


% --- Executes during object creation, after setting all properties.
function PositionX0_CreateFcn(hObject, eventdata, handles) %#ok

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in blankstim.
function blankstim_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to blankstim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of blankstim


function PositionY0_Callback(hObject, eventdata, handles) %#ok
StimType_Callback(handles.StimType,eventdata,handles);

% --- Executes during object creation, after setting all properties.
function PositionY0_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to PositionY0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in FullFlicker.
function FullFlicker_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to FullFlicker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FullFlicker


function phasePeriod_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to phasePeriod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of phasePeriod as text
%        str2double(get(hObject,'String')) returns contents of phasePeriod as a double


% --- Executes during object creation, after setting all properties.
function phasePeriod_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to phasePeriod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function nReps_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to nReps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nReps as text
%        str2double(get(hObject,'String')) returns contents of nReps as a double


% --- Executes during object creation, after setting all properties.
function nReps_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to nReps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function stimulusGroups_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to stimulusGroups (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stimulusGroups as text
%        str2double(get(hObject,'String')) returns contents of stimulusGroups as a double

% --- Executes during object creation, after setting all properties.
function stimulusGroups_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to stimulusGroups (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in squaregratings.
function squaregratings_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to squaregratings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of squaregratings


function updateparamsfrommovie(fullmoviename,handles)

% update variables after select
load(fullmoviename,'MovieMag','MovieRate','period_sec');

if exist('MovieMag','var')
    % movie meta-data variables saved:
    %  'duration_sec','period_sec','MovieMag','MovieRate','screenDistanceCm'

    set(handles.MovieMag,'String',num2str(MovieMag)); %#ok
    set(handles.MovieRate,'String',num2str(MovieRate)); %#ok
    set(handles.phasePeriod,'String',num2str(period_sec)); %#ok

    % other variables that could be saved:
    % duration_sec (default is to play whole movie)
    % screenDistanceCm (irrelevant for movies)
    % xsize (no ysize)
end



% --- Executes on button press in eyeCond0.
function eyeCond0_Callback(hObject, eventdata, handles)
% hObject    handle to eyeCond0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function eyeCond0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eyeCond0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function maskcenterx_Callback(hObject, eventdata, handles)
% hObject    handle to maskcenterx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maskcenterx as text
%        str2double(get(hObject,'String')) returns contents of maskcenterx as a double


% --- Executes during object creation, after setting all properties.
function maskcenterx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maskcenterx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maskcentery_Callback(hObject, eventdata, handles)
% hObject    handle to maskcentery (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maskcentery as text
%        str2double(get(hObject,'String')) returns contents of maskcentery as a double


% --- Executes during object creation, after setting all properties.
function maskcentery_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maskcentery (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maskradiusx_Callback(hObject, eventdata, handles)
% hObject    handle to maskradiusx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maskradiusx as text
%        str2double(get(hObject,'String')) returns contents of maskradiusx as a double


% --- Executes during object creation, after setting all properties.
function maskradiusx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maskradiusx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maskradiusy_Callback(hObject, eventdata, handles)
% hObject    handle to maskradiusy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maskradiusy as text
%        str2double(get(hObject,'String')) returns contents of maskradiusy as a double


% --- Executes during object creation, after setting all properties.
function maskradiusy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maskradiusy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maskcenterdeg_Callback(hObject, eventdata, handles)
% hObject    handle to maskcenterdeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maskcenterdeg as text
%        str2double(get(hObject,'String')) returns contents of maskcenterdeg as a double


% --- Executes during object creation, after setting all properties.
function maskcenterdeg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maskcenterdeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maskmeanradiusdeg_Callback(hObject, eventdata, handles)
% hObject    handle to maskmeanradiusdeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maskmeanradiusdeg as text
%        str2double(get(hObject,'String')) returns contents of maskmeanradiusdeg as a double


% --- Executes during object creation, after setting all properties.
function maskmeanradiusdeg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maskmeanradiusdeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in popmenuMask.
function popmenuMask_Callback(hObject, eventdata, handles)
% hObject    handle to popmenuMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popmenuMask contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popmenuMask


% --- Executes during object creation, after setting all properties.
function popmenuMask_CreateFcn(hObject, eventdata,handles)
% hObject    handle to popmenuMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function masktex =  helperMakeMask(rx,ry,window,masktype,cnow,handles)
% BA
tic
[width, height]=Screen('WindowSize', window);

width = width*2.1;% height a width multipled by 2.1 to make mask can be moved anywhere in the
height = height*2.1;
% screen and still mask
white = WhiteIndex(window);
black = BlackIndex(window);
grey = round(0.5*(black+white));

% We create a Luminance+Alpha matrix for use as transparency mask:
[x,y]=meshgrid([1:width]-width/2,[1:height]-height/2);
% Layer 1 (Luminance) is filled with luminance value 'gray' of the
% background.
maskimg=ones(height,width,2) * grey;
% Layer 2 (Transparency aka Alpha) is filled with gaussian transparency
% mask.

switch masktype % note param.nMask should equal to the number of cases
    case 0 % gaussianm
        maskimg(:,:,2)=255 - exp(-((x/rx).^2)-((y/ry).^2))*255;
    case 1 % eliptical aperature
        maskimg(:,:,2) = 255;
        maskimg((height/2-rx):(height/2+rx-1),(width/2-ry):(width/2+ry-1),2)= (~makeElipse(rx,ry))*255;
    case 2 % inverted eliptical aperature
        maskimg(:,:,2) = 0;
        maskimg((height/2-rx):(height/2+rx-1),(width/2-ry):(width/2+ry-1),2)= (makeElipse(rx,ry))*255;
    case 3 % no mask
        maskimg(:,:,2) = 0;
%     case 4 % Working on cross grating
%         Duration = str2double(get(handles.Duration,'String'));
%         FrameHz = round(str2double(get(handles.FrameHz,'String')));
% 
%         cnow = UCparams.c;
% 
%         [frm]= generateGratings_blit(handles.orient(cnow),handles.freq(cnow),handles.TempFreq(cnow),handles.phase(cnow),handles.contrast(cnow),1/FrameHz,  handles.degPerPix,width,height,FrameHz,black,white);
%         maskimg(:,:,2) = frm;
%         maskimg((height/2-rx):(height/2+rx-1),(width/2-ry):(width/2+ry-1),2)= (~makeElipse(rx,ry))*255;
% 
end
% Build a single transparency mask texture:
masktex=Screen('MakeTexture', window, maskimg);
% masktex=Screen('MakeTexture', window, squeeze(frm));

% Screen('DrawTexture', window,masktex)
% Screen('Flip',window);
toc


function UCkeys = declareUCkeys()
% Set keys.
% UCkeys.rightKey = KbName('RightArrow');
% UCkeys.leftKey = KbName('LeftArrow');
UCkeys.commaKey = KbName('.>');
UCkeys.periodKey = KbName(',');
UCkeys.upKey = KbName('UpArrow');
UCkeys.downKey = KbName('DownArrow');
UCkeys.mKey = KbName('m');
UCkeys.spaceKey = KbName('space');
UCkeys.aKey = KbName('a');  % auto bautoChangeContrast
UCkeys.lKey = KbName('l');
UCkeys.escKey = KbName('ESCAPE');
UCkeys.buttons = 0; % When the user clicks the mouse, 'buttons' becomes nonzero.
UCkeys.strikeKey = KbName('`');

function UCparams = helperUserControl(keyCode,UCkeys,UCparams,params,handles)
% function that handles key strokes during stimulus presentation
% UCkeys - user controlled keys
% UCparams - user controlled params (may be changed in this funciton)
% params - should not be changed in this function
%
% BA

% ADD rotation of texture

if keyCode(UCkeys.spaceKey) % next condition
    UCparams.c = UCparams.c+1;
    if UCparams.c > params.nCond
        UCparams.c = 1;
    end
elseif keyCode(UCkeys.commaKey) % rotate bar
    UCparams.rotation = mod(UCparams.rotation + 15,360);
elseif keyCode(UCkeys.periodKey) % rotate bar
    UCparams.rotation = mod(UCparams.rotation - 15,360);
elseif keyCode(UCkeys.upKey)
    UCparams.rx = min(UCparams.rx+params.rStep,params.RMAX);
    UCparams.ry = UCparams.rx; % ADD so that it doesn't have to creat new mask
    UCparams.masktex =  helperMakeMask(UCparams.rx,UCparams.ry,params.window,UCparams.masktype,UCparams.c,handles);
elseif keyCode(UCkeys.downKey)
    UCparams.rx = max(UCparams.rx-params.rStep,0);
    UCparams.ry = UCparams.rx;
    UCparams.masktex =  helperMakeMask(UCparams.rx,UCparams.ry,params.window,UCparams.masktype,UCparams.c,handles);
elseif keyCode(UCkeys.mKey)
    UCparams.masktype = UCparams.masktype+1;
    UCparams.masktype = mod(UCparams.masktype,params.nMasks);
    % to DO replace so that mask doesn't have to be
    % recalcualted
    UCparams.masktex =  helperMakeMask(UCparams.rx,UCparams.ry,params.window,UCparams.masktype,UCparams.c,handles);
elseif  keyCode(UCkeys.lKey) % toggle lock mask location
    UCparams.lockMask = ~UCparams.lockMask;
    if UCparams.lockMask
        showCursor;
    else
        hideCursor;
    end
elseif  keyCode(UCkeys.aKey) % toggle lock mask location
    UCparams.bautoChangeContrast = ~UCparams.bautoChangeContrast;
elseif  keyCode(UCkeys.escKey) % end after this condition finishes
    UCparams.doneStim = 1;
    keyspressed = KbName(find(keyCode));
    disp(sprintf('Exit on %s key pressed',keyspressed));
elseif  keyCode(UCkeys.strikeKey) % end immediately
    UCparams.break= 1;
    UCparams.doneStim = 1;
    keyspressed = KbName(find(keyCode));
    disp(sprintf('Exit on %s key pressed',keyspressed));
end

