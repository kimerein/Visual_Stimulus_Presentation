try
    AssertOpenGL;
    screenNumber=0; 
    window=Screen('OpenWindow',screenNumber,0,[],[],2);
    stepsize=10;
    steps=ceil(255/stepsize);
    for i=0:steps
        pixelvalue=i*stepsize;
        if pixelvalue>255
            pixelvalue=255;
        end
        Screen('FillRect',window, pixelvalue);
        Screen('DrawText',window,num2str(pixelvalue),0,0,255);
        Screen('Flip', window);
        [x,y,buttons] = GetMouse;
        while any(buttons) % wait for release
            [x,y,buttons] = GetMouse;
        end
        while 1
            [mX, mY, buttons] = GetMouse;
            if buttons(1)
                break;
            end
        end   
    end
    Screen('CloseAll');
catch
    Screen('CloseAll');
    psychrethrow(psychlasterror);
end