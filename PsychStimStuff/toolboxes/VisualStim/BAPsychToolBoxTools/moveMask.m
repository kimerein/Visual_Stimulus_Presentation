function moveMask% display a bar in psychtoolbox
try
    screenNumber = 1;

    KbName('UnifyKeyNames');
    FrameHz=Screen('FrameRate',screenNumber);


    [window,windowRect]=Screen(screenNumber,'OpenWindow',0);   %%% open grey window
    % Enable alpha blending with proper blend-function. We need it
    % for drawing of our alpha-mask (gaussian aperture):
    Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


    white = WhiteIndex(window);
    black = BlackIndex(window);
    grey = round(0.5*(black+white));
    %
    % 	Screen('FillRect',window, gray);
    % 	Screen('Flip', window);
    % 	Screen('FillRect',window, gray);


    imageRect = windowRect;
    Screen('DrawText',window,sprintf('Generating stimuli'),10,30,white);
    Screen('Flip',window);

  freq = [0.5 0.15 0.08 0.025];
  
    orient = [0 0 0 0];
    speed = [25 25 25 25];
    contrast = [1 1 1 1];
    length = [0 0 0 0];
    positionX =[-5 -5 -5 -5];
    Duration = 1;
    degPerPix = .0716; % REPLACe
    sizeLut = 256;

    % grating specific params
    TempFreq =1;
    phase0 = 0;
    Rot = 90;
    
    % mask radii
    rx = 90; % in pixels %% add change to degrs
    ry = 90;
    rStep = 20;
    RMAX = max(imageRect);
    lockMask = 1;

    masktype = 1;
    masktex =  helperMakeMask(rx,ry,window,masktype);

    
    
    buttons = 0; % When the user clicks the mouse, 'buttons' becomes nonzero.
    mX = imageRect(3)/2; % The x-coordinate of the mouse cursor
    mY = imageRect(4)/2; % The y-coordinate of the mouse cursor

    nCond = size(orient,2);

    % make textures
    textures = zeros(nCond,1);
    for c = 1:nCond
        % change to texture (so can rotate)
        %     [img cl] = generateBars_lut(orient(c),freq(c),speed(c),contrast(c),length(c), positionX(c),Duration, degPerPix,imageRect(3),imageRect(4),FrameHz,black,white,sizeLut);
        [img cl] = generateGratings_lut(orient(c),freq(c),TempFreq,phase0,contrast(c),Duration, degPerPix,imageRect(3),imageRect(4),frameRate, black, white,sizeLut);
        
        fprintf('done generating')
        if c==1
            cluts = zeros(256,3,size(cl,3),nCond);
        end
        cluts(:,:,:,c)=floor(cl);
        textures(c,1)=Screen('MakeTexture',window,img);
    end
    FrameWait=1; % BA?


    ifi_duration = Screen('GetFlipInterval', window);

    % Set keys.
    rightKey = KbName('RightArrow');
    leftKey = KbName('LeftArrow');
    upKey = KbName('UpArrow');
    downKey = KbName('DownArrow');
    mKey = KbName('m');
    cKey = KbName('c');
    lKey = KbName('l');
    escKey = KbName('ESCAPE');
    
    %     escapeKey = KbName('ESCAPE');

    doneStim = 0;
    texturec = 1;
    lastsecs = [];

    c=1;
            vbl = Screen('Flip',window);  %%% initial flip, to sync with vertical blank

    while ~doneStim && ~any(buttons)
        priorityLevel = MaxPriority(window);
        Priority(priorityLevel);
        %     clut;

        clutcond = (squeeze(cluts(:,:,:,c)));
        nFrames =  size(cluts,3);
        
        for f = 1:nFrames
            if ~lockMask
                [mX, mY, buttons] = GetMouse;
            end
            moglClutBlit(window,textures(c),clutcond(:,:,f));
            Screen('DrawTexture', window, masktex, [],CenterRectOnPoint(windowRect*2, mX, mY));
            %% if needed ADD rotation here

            vbl=Screen('Flip', window, vbl + (1 - 0.5) * ifi_duration);
            


            % look for stop message on keyboard
            
            [keyIsDown, secs, keyCode] = KbCheck;
            if isempty(lastsecs)|( secs -lastsecs) >= 1/10; % set a maximum rate at which a key can be pressed otherwise holding button makes bar spin to fast
                if keyIsDown %%% charavail would be much better, but doesn't seem to work
                    lastsecs= secs;
                    if keyCode(cKey)
                        c = c+1;
                        if c > nCond
                            c = 1;
                        end
                    elseif keyCode(upKey)
                        rx = min(rx+rStep,RMAX);
                        ry = rx; % ADD so that it doesn't have to creat new mask
                        masktex =  helperMakeMask(rx,ry,window,masktype);
                    elseif keyCode(downKey)
                        rx = max(rx-rStep,0);
                        ry = rx;
                        masktex =  helperMakeMask(rx,ry,window,masktype);
                    elseif keyCode(mKey)
                        masktype = masktype+1;
                        masktype = mod(masktype,4);
                        % to DO replace so that mask doesn't have to be
                        % recalcualted
                        masktex =  helperMakeMask(rx,ry,window,masktype);
                    elseif  keyCode(lKey) % toggle lock mask location
                        lockMask = ~lockMask;
                    elseif  keyCode(escKey)
                        doneStim = 1;
                        keyspressed = KbName(find(keyCode));
                        disp(sprintf('Exit on %s key pressed',keyspressed{1}));
                    end
                end
            

            end
        end
    end


    % clean up
    moglClutBlit;
    ListenChar(1);
    %      Screen('LoadNormalizedGammaTable',window,flat_clut);
    screen('CloseAll');

catch
    ShowCursor;
    moglClutBlit;
    ListenChar(1);
    %      Screen('LoadNormalizedGammaTable',window,flat_clut);
    screen('CloseAll');
end

function masktex =  helperMakeMask(rx,ry,window,masktype)
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

switch masktype
    case 0 % gaussian
        maskimg(:,:,2)=255 - exp(-((x/rx).^2)-((y/ry).^2))*255;
    case 1 % eliptical aperature
        maskimg(:,:,2) = 255;
        maskimg((height/2-rx):(height/2+rx-1),(width/2-ry):(width/2+ry-1),2)= (~makeElipse(rx,ry))*255;
    case 2 % inverted eliptical aperature
        maskimg(:,:,2) = 0;
        maskimg((height/2-rx):(height/2+rx-1),(width/2-ry):(width/2+ry-1),2)= (makeElipse(rx,ry))*255;
    case 3 % no mask
        maskimg(:,:,2) = 0;
end
% Build a single transparency mask texture:
masktex=Screen('MakeTexture', window, maskimg);
toc

