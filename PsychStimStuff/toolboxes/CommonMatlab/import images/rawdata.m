
%   Data File Conversion Program for converting Multiphoton data files
%   Converts MPD files and older Labview Dat files to multiframe .tif
%       and/or .mat files under the variable name 'mydata'
%
%Written by Phil Tsai 02/14/04
%Last Revision : 02/26/08

%Choose File to Read
[fileName,pathName] = uigetfile('*.dat; *.mpd','Choose a data file to convert');
inputFileName = [pathName,fileName];
fileNameLength = length(inputFileName);
fileRoot = inputFileName(1,1:fileNameLength-4);
fileExt = inputFileName(fileNameLength-2:fileNameLength);
cd(pathName)

fileType = '';
if strcmpi(fileExt,'mpd'), fileType = '*.mpd'; end
if strcmpi(fileExt,'dat'), fileType = '*.dat'; end


%Checkbox Dialog Box for Processing Options
keepGoingCheckbox = 0;
button1Val = 1; %Save to 16-bit multiframe tiff
button2Val = 0; %Tile channels together
button3Val = 0; %Save to .mat file as mydata
button4Val = 1; %Convert all files in folder
button5Val = 1; %Skip previously converted files
button6Val = 0; %Show waitbar
%-----
while keepGoingCheckbox==0; 
  hitme = 0;
  f55 = figure(55); clf(55);
  windowHeight = 190;
  set(f55,'Name','Processing Options',...
     'MenuBar','none','NumberTitle','off',...
     'Position',[350,200,240,windowHeight],...
     'Color',[0.2,0.8,0.2]);
  button1 = uicontrol(f55,'Style','checkbox',...
     'BackgroundColor',[0.2,0.8,0.2],...
     'ForegroundColor',[0,0,0],...
     'Value',button1Val==1,...
     'Position',[25,windowHeight-30,200,15],...
     'String','Save to 16-bit multiframe tiff',...
     'Callback','hitme=1;button1Val=~button1Val;');
  button2 = uicontrol(f55,'Style','checkbox',...
     'BackgroundColor',[0.2,0.8,0.2],...
     'ForegroundColor',[0,0,0],...
     'Value',button2Val==1,...
     'Position',[25,windowHeight-50,200,15],...
     'String','Tile Channels Together',...
     'Callback','hitme=1;button2Val=~button2Val;'); 
  button3 = uicontrol(f55,'Style','checkbox',...
     'BackgroundColor',[0.2,0.8,0.2],...
     'ForegroundColor',[0,0,0],...
     'Value',button3Val==1,...
     'Position',[25,windowHeight-70,200,15],...
     'String','Save as .mat file (variable = mydata)',...
     'Callback','hitme=1;button3Val=~button3Val;'); 
  button4 = uicontrol(f55,'Style','checkbox',...
     'BackgroundColor',[0.2,0.8,0.2],...
     'ForegroundColor',[0,0,0],...
     'Value',button4Val==1,...
     'Position',[25,windowHeight-90,200,15],...
     'String','Convert all files in folder',...
     'Callback','hitme=1;button4Val=~button4Val;'); 
   button5 = uicontrol(f55,'Style','checkbox',...
     'BackgroundColor',[0.2,0.8,0.2],...
     'ForegroundColor',[0,0,0],...
     'Value',button5Val==1,...
     'Position',[25,windowHeight-110,200,15],...
     'String','Skip previously-converted files',...
     'Callback','hitme=1;button5Val=~button5Val;'); 
   button6 = uicontrol(f55,'Style','checkbox',...
     'BackgroundColor',[0.2,0.8,0.2],...
     'ForegroundColor',[0,0,0],...
     'Value',button6Val==1,...
     'Position',[25,windowHeight-130,200,15],...
     'String','Show Waitbar',...
     'Callback','hitme=1;button6Val=~button6Val;'); 
  finishButton = uicontrol(f55,'Style','pushbutton',...
     'BackgroundColor',[0.2,0.6,0.2],...
     'Value',keepGoingCheckbox,...
     'Position',[100,windowHeight-175,50,25],...
     'String','OK',...
     'Callback','hitme=1; keepGoingCheckbox=1;');
  while hitme == 0;
     pause(0.1);
  end; %while hitme
end %while keepGoingCheckbox 
%-----
convertToTiffFlag = button1Val;
tileChannelsFlag = button2Val;
storeMatFileFlag = button3Val;
convertAllFlag = button4Val;
skipConvertedFlag = button5Val;
waitbarFlag = button6Val;
close(f55);
%----------
%If option is selected, get all files in directory
if convertAllFlag == 1,
    filesToConvert = dir(fileType);
else %if convertAllFlag == 1,
    filesToConvert.name = fileName;
end %if convertAllFlag == 1,
numFiles = length(filesToConvert);
%----------
%Now begin the appropriate extraction alogrithim
switch fileType
    case '*.mpd' %MPScan .mpd file
        %---------
        f287 = figure(287);
        close(f287);
        f287 = figure(287);
        set(f287,'Visible','Off')
        %----------
        for fileIter = 1:numFiles
            currentFileName = filesToConvert(fileIter).name;
            fileRoot = currentFileName(1,1:length(currentFileName)-4);
            progressReport = [' File #', num2str(fileIter)];
            progressReport = [progressReport, ' of ',num2str(numFiles)];
            %---
            mpfile = actxcontrol('MPfile.Data',[0,0,500,500],f287);
            OpenFile = invoke(mpfile,'OpenMPFile',currentFileName);
            %---
            %Read all the header information, but currently, do nothing with most of it
            Header.Scan_Mode = invoke(mpfile,'ReadParameter','Scan Mode');
            Header.Frame_Width = invoke(mpfile,'ReadParameter','Frame Width');
            Header.Frame_Height = invoke(mpfile,'ReadParameter','Frame Height');
            Header.Frame_Count = invoke(mpfile,'ReadParameter','Frame Count');
            Header.X_Position = invoke(mpfile,'ReadParameter','X Position');
            Header.Y_Position = invoke(mpfile,'ReadParameter','Y Position');
            Header.Z_Position = invoke(mpfile,'ReadParameter','Z Position');
            Header.Stack_Count = invoke(mpfile,'ReadParameter','Stack Count');
            Header.Z_Interval = invoke(mpfile,'ReadParameter','z- Interval');
            Header.Averaging_Count = invoke(mpfile,'ReadParameter','Averaging Count');
            Header.Repeat_Count = invoke(mpfile,'ReadParameter','Repeat Count');
            Header.Magnification = invoke(mpfile,'ReadParameter','Magnification');
            Header.Rotation = invoke(mpfile,'ReadParameter','Rotation');
            Header.X_Frame_Offset = invoke(mpfile,'ReadParameter','X Frame Offset');
            Header.Y_Frame_Offset = invoke(mpfile,'ReadParameter','Y Frame Offset');
            Header.Channel_Name1 = invoke(mpfile,'ReadParameter','Channel Name (1)');
            Header.Channel_Name2 = invoke(mpfile,'ReadParameter','Channel Name (2)');
            Header.Channel_Name3 = invoke(mpfile,'ReadParameter','Channel Name (3)');
            Header.Channel_Name4 = invoke(mpfile,'ReadParameter','Channel Name (4)');
            Header.Enabled1 = invoke(mpfile,'ReadParameter','Enabled (1)');
            Header.Enabled2 = invoke(mpfile,'ReadParameter','Enabled (2)');
            Header.Enabled3 = invoke(mpfile,'ReadParameter','Enabled (3)');
            Header.Enabled4 = invoke(mpfile,'ReadParameter','Enabled (4)');
            Header.Input_Range1 = invoke(mpfile,'ReadParameter','Input Range (1)');
            Header.Input_Range2 = invoke(mpfile,'ReadParameter','Input Range (2)');
            Header.Input_Range3 = invoke(mpfile,'ReadParameter','Input Range (3)');
            Header.Input_Range4 = invoke(mpfile,'ReadParameter','Input Range (4)');
            Header.Channel_Unit3 = invoke(mpfile,'ReadParameter','Channel Unit (3)');
            Header.Channel_Unit4 = invoke(mpfile,'ReadParameter','Channel Unit (4)');
            Header.Channel_Prefix3 = invoke(mpfile,'ReadParameter','Channel Prefix (3)');
            Header.Channel_Prefix4 = invoke(mpfile,'ReadParameter','Channel Prefix (4)');
            Header.Conversion_Factor3 = invoke(mpfile,'ReadParameter','Conversion Factor (3)');
            Header.Conversion_Factor4 = invoke(mpfile,'ReadParameter','Conversion Factor (4)');
            Header.Offset3 = invoke(mpfile,'ReadParameter','Offset (3)');
            Header.Offset4 = invoke(mpfile,'ReadParameter','Offset (4)');
            Header.Data_Point_Per_Frame3 = invoke(mpfile,'ReadParameter','Data Point Per Frame (3)');
            Header.Data_Point_Per_Frame4 = invoke(mpfile,'ReadParameter','Data Point Per Frame (4)');
            Header.Comments = invoke(mpfile,'ReadParameter','Comments');
            %---
            numFrames = str2num(Header.Frame_Count);
            xsize = str2num(Header.Frame_Width);
            ysize = str2num(Header.Frame_Height);
            Ch1Flag = 0;
            Ch2Flag = 0;
            Ch3Flag = 0;
            Ch4Flag = 0;
            if strcmp(Header.Enabled1,'True'), Ch1Flag = 1; numChannels = 1; end
            if strcmp(Header.Enabled2,'True'), Ch2Flag = 1; numChannels = 2; end
            if strcmp(Header.Enabled3,'True'), Ch3Flag = 1; numChannels = 3; end
            if strcmp(Header.Enabled4,'True'), Ch4Flag = 1; numChannels = 4; end
            %---
            if storeMatFileFlag == 1,
                mydata = zeros([ysize,xsize,numFrames],'uint16');
            end%if storeMatFileFlag == 1,
            %---
            if Ch1Flag == 1,
                outputTifName1 = [fileRoot,'-Ch1.tif'];
                outputMatName1 = [fileRoot,'-Ch1.mat'];
                tiledTifName = [fileRoot,'-ChTiled.tif'];
                tifExistsFlag = exist(outputTifName1,'file');
                matExistsFlag = exist(outputMatName1,'file');
                tiledExistsFlag = exist(tiledTifName,'file');
                fileExistsFlag = 1;
                if and(convertToTiffFlag==1 , tifExistsFlag==0), fileExistsFlag = 0; end
                if and(storeMatFileFlag==1 , matExistsFlag==0), fileExistsFlag = 0; end  
                if and(tileChannelsFlag==1 , tiledExistsFlag==0), fileExistsFlag = 0; end
                if or(skipConvertedFlag==0, fileExistsFlag==0),
                    if waitbarFlag == 1, wb = waitbar(0,['Converting Ch.1, ',progressReport]); end
                    writeMode = 'overwrite';
                    for currentFrame = 1:numFrames
                        dataFrame = mpfile.ReadFrameData(1,currentFrame);
                        imageFrame = transpose(reshape(dataFrame,xsize,ysize));
                        if and(convertToTiffFlag == 1 , tifExistsFlag == 0);
                            imwrite(uint16(imageFrame),outputTifName1,'tiff','Compression','none','WriteMode',writeMode,'Resolution',1);
                            writeMode = 'append';
                        end%if and(convertToTiffFlag == 1 , tifExistsFlag == 0);
                        %---
                        if or(storeMatFileFlag == 1, tileChannelsFlag == 1),
                            mydata(:,:,currentFrame) = imageFrame;
                        end%if or(storeMatFileFlag == 1, tileChannelsFlag == 1),
                        %---
                        if waitbarFlag == 1, waitbar(currentFrame/numFrames,wb); end
                    end %for currentFrame
                    if waitbarFlag == 1,close(wb);end
                    %---
                    if storeMatFileFlag == 1,
                        save(outputMatName1,'mydata');
                    end%if storeMatFileFlag == 1,
                    %---
                    if tileChannelsFlag == 1,
                        dataCh1 = mydata;
                    end%if tileChannelsFlag == 1,
                end%%if or(skipConvertedFlag==0, fileExistsFlag==0),
            end %if Ch1Flag == 1;
            %-----
            if Ch2Flag == 1,
                outputTifName2 = [fileRoot,'-Ch2.tif'];
                outputMatName2 = [fileRoot,'-Ch2.mat'];
                tiledTifName = [fileRoot,'-ChTiled.tif'];
                tifExistsFlag = exist(outputTifName2,'file');
                matExistsFlag = exist(outputMatName2,'file');
                tiledExistsFlag = exist(tiledTifName,'file');
                fileExistsFlag = 1;
                if and(convertToTiffFlag==1 , tifExistsFlag==0), fileExistsFlag = 0; end
                if and(storeMatFileFlag==1 , matExistsFlag==0), fileExistsFlag = 0; end  
                if and(tileChannelsFlag==1 , tiledExistsFlag==0), fileExistsFlag = 0; end
                if or(skipConvertedFlag==0, fileExistsFlag==0),
                    if waitbarFlag == 1, wb = waitbar(0,['Converting Ch.2, ',progressReport]); end
                    writeMode = 'overwrite';
                    for currentFrame = 1:numFrames
                        dataFrame = mpfile.ReadFrameData(2,currentFrame);
                        imageFrame = transpose(reshape(dataFrame,xsize,ysize));
                        if and(convertToTiffFlag == 1 , tifExistsFlag == 0);
                            imwrite(uint16(imageFrame),outputTifName2,'tiff','Compression','none','WriteMode',writeMode,'Resolution',1);
                            writeMode = 'append';
                        end%if and(convertToTiffFlag == 1 , tifExistsFlag == 0);
                        %---
                        if or(storeMatFileFlag == 1, tileChannelsFlag == 1),
                            mydata(:,:,currentFrame) = imageFrame;
                        end%if or(storeMatFileFlag == 1, tileChannelsFlag == 1),
                        %---
                        if waitbarFlag == 1, waitbar(currentFrame/numFrames,wb); end
                    end %for currentFrame
                    if waitbarFlag == 1,close(wb);end
                    %---
                    if storeMatFileFlag == 1,
                        save(outputMatName2,'mydata');
                    end%if storeMatFileFlag == 1,
                    %---
                    if tileChannelsFlag == 1,
                        dataCh2 = mydata;
                    end%if tileChannelsFlag == 1,
                end%if or(skipConvertedFlag==0, fileExistsFlag==0),
            end %if Ch2Flag == 1;
            %-----
            if Ch3Flag == 1,
                outputTifName3 = [fileRoot,'-Ch3.tif'];
                outputMatName3 = [fileRoot,'-Ch3.mat'];
                tiledTifName = [fileRoot,'-ChTiled.tif'];
                tifExistsFlag = exist(outputTifName3,'file');
                matExistsFlag = exist(outputMatName3,'file');
                tiledExistsFlag = exist(tiledTifName,'file');
                fileExistsFlag = 1;
                if and(convertToTiffFlag==1 , tifExistsFlag==0), fileExistsFlag = 0; end
                if and(storeMatFileFlag==1 , matExistsFlag==0), fileExistsFlag = 0; end  
                if and(tileChannelsFlag==1 , tiledExistsFlag==0), fileExistsFlag = 0; end
                if or(skipConvertedFlag==0, fileExistsFlag==0),
                    if waitbarFlag == 1, wb = waitbar(0,['Converting Ch.3, ',progressReport]);end
                    writeMode = 'overwrite';
                    for currentFrame = 1:numFrames
                        dataFrame = mpfile.ReadFrameData(3,currentFrame);
                        imageFrame = transpose(reshape(dataFrame,xsize,ysize));
                        if and(convertToTiffFlag == 1 , tifExistsFlag == 0);
                            imwrite(uint16(imageFrame),outputTifName3,'tiff','Compression','none','WriteMode',writeMode,'Resolution',1);
                            writeMode = 'append';
                        end%if and(convertToTiffFlag == 1 , tifExistsFlag == 0);
                        %---
                        if or(storeMatFileFlag == 1, tileChannelsFlag == 1),
                            mydata(:,:,currentFrame) = imageFrame;
                        end%if or(storeMatFileFlag == 1, tileChannelsFlag == 1),
                        %---
                        if waitbarFlag == 1, waitbar(currentFrame/numFrames,wb); end
                    end %for currentFrame
                    if waitbarFlag == 1,close(wb);end
                    %---
                    if storeMatFileFlag == 1,
                        save(outputMatName3,'mydata');
                    end%if storeMatFileFlag == 1,
                    %---
                    if tileChannelsFlag == 1,
                        dataCh3 = mydata;
                    end%if tileChannelsFlag == 1,
                end%if or(skipConvertedFlag==0, fileExistsFlag==0),
            end %if Ch3Flag == 1;
            %-----
            if Ch4Flag == 1,
                outputTifName4 = [fileRoot,'-Ch4.tif'];
                outputMatName4 = [fileRoot,'-Ch4.mat'];
                tiledTifName = [fileRoot,'-ChTiled.tif'];
                tifExistsFlag = exist(outputTifName4,'file');
                matExistsFlag = exist(outputMatName4,'file');
                tiledExistsFlag = exist(tiledTifName,'file');
                fileExistsFlag = 1;
                if and(convertToTiffFlag==1 , tifExistsFlag==0), fileExistsFlag = 0; end
                if and(storeMatFileFlag==1 , matExistsFlag==0), fileExistsFlag = 0; end  
                if and(tileChannelsFlag==1 , tiledExistsFlag==0), fileExistsFlag = 0; end
                if or(skipConvertedFlag==0, fileExistsFlag==0),
                    if waitbarFlag == 1, wb = waitbar(0,['Converting Ch.4, ',progressReport]); end
                    writeMode = 'overwrite';
                    for currentFrame = 1:numFrames
                        dataFrame = mpfile.ReadFrameData(4,currentFrame);
                        imageFrame = transpose(reshape(dataFrame,xsize,ysize));
                        if and(convertToTiffFlag == 1 , tifExistsFlag == 0);
                            imwrite(uint16(imageFrame),outputTifName4,'tiff','Compression','none','WriteMode',writeMode,'Resolution',1);
                            writeMode = 'append';
                        end%if and(convertToTiffFlag == 1 , tifExistsFlag == 0);
                        %---
                        if or(storeMatFileFlag == 1, tileChannelsFlag == 1),
                            mydata(:,:,currentFrame) = imageFrame;
                        end%if or(storeMatFileFlag == 1, tileChannelsFlag == 1),
                        %---
                        if waitbarFlag == 1, waitbar(currentFrame/numFrames,wb); end
                    end %for currentFrame
                    if waitbarFlag == 1,close(wb);end
                    %---
                    if storeMatFileFlag == 1,
                        save(outputMatName4,'mydata');
                    end%if storeMatFileFlag == 1,
                    %---
                    if tileChannelsFlag == 1,
                        dataCh4 = mydata;
                    end%if tileChannelsFlag == 1,
                end%if or(skipConvertedFlag==0, fileExistsFlag==0),
            end %if Ch4Flag == 1;
            
            %TileChannelsTogether if desired
            if tileChannelsFlag == 1,
                tiledTifName = [fileRoot,'-ChTiled.tif'];
                fileExistsFlag = exist(tiledTifName,'file');
                if or(skipConvertedFlag==0, fileExistsFlag==0),
                    tiledStackSize = [ysize,numChannels*xsize,numFrames];
                    tiledStack = zeros(tiledStackSize,'uint16');
                    if Ch1Flag==1,
                        tiledStack(:,1:xsize,:) = dataCh1;
                    end%if Ch1Flag==1,
                    if Ch2Flag==1,
                        tiledStack(:,xsize+1:2*xsize,:) = dataCh2;
                    end%if Ch2Flag==1,
                    if Ch3Flag==1,
                        tiledStack(:,2*xsize+1:3*xsize,:) = dataCh3;
                    end%if Ch3Flag==1,
                    if Ch4Flag==1,
                        tiledStack(:,3*xsize+1:4*xsize,:) = dataCh4;
                    end%if Ch4Flag==1,
                    if waitbarFlag == 1, wb = waitbar(0,['Tiling Channels, ',progressReport]);end
                    writeMode = 'overwrite';
                    for currentFrame = 1:numFrames
                        imwrite(uint16(tiledStack(:,:,currentFrame)),tiledTifName,'tiff','Compression','none','WriteMode',writeMode,'Resolution',1);
                        writeMode = 'append';
                        if waitbarFlag == 1, waitbar(currentFrame/numFrames,wb); end
                    end%for currentFrame = 1:numFrames
                    if waitbarFlag == 1,close(wb);end
                    clear tiledStack
                end%if or(skipConvertedFlag==0, fileExistsFlag==0),
                clear dataCh*
            end%if tileChannelsFlag == 1,
        end%for fileIter = 1:numfiles
        close(f287); 
    %----------
    
    case '*.dat' %Labview .dat file
        for fileIter = 1:numFiles,
            currentFileName = filesToConvert(fileIter).name;
            fileRoot = currentFileName(1,1:length(currentFileName)-4);
            
            %Read the headers into a structure with the grabframeheaders function file
            [myheaders,mychannels,num_datframes,header_sizes] = grabframeheaders(['\',currentFileName]);%Open the file

            %Determine # of frames at beginning of each channel with same frame
            %size using the samesizeframes function file
            samesize = samesizeframes(myheaders);


            %Start by figuring (from the headers) the number and size of the
            % frames in each channel to be processed
            whichChannel = find(samesize>0);
            numPics = min(samesize(whichChannel)); %num of pics to process in per channel
            numChannels = sum(samesize>0);
            numFrames = numPics*numChannels; %total num of frames to analyze


            %Go through all channels and create an ouput file of each selected type for
            % each channel; write the appropriate headers. 
            for currentChannelPointer = 1:numChannels, %currentChannelPointer = nth loop
                currentChannel = whichChannel(currentChannelPointer); %current channel = actual channel ID
                %get frame parameters for each channel from first example of each channel
                gain = str2num(myheaders(currentChannel+1).gain);
                gain_offset = str2num(myheaders(currentChannel+1).gain_offset);
                xsize = str2num(myheaders(currentChannel+1).xsize)-str2num(myheaders(currentChannel+1).offset);
                ysize = str2num(myheaders(currentChannel+1).ysize);
                frame_size(currentChannel) = xsize*ysize;
                if samesize(currentChannel)~=0, %if there is data in this channel create an output file
                    if convertToTiffFlag==1,
                        outputTiffFileName{currentChannelPointer}=[fileRoot,'-Ch',num2str(currentChannel),'.tif'];
                    end%if convertToTiffFlag==1
                end %if samesize(currentChannel)~=0
            end %for currentChannelPointer = 1:numChannels

            %Start extracting the images in sequential order
            %and writing them to files frame by frame
            fidr = fopen(currentFileName,'r','ieee-be');
            channel_stub = find(samesize~=0);
            channel_seq = repmat(channel_stub,[1,(max(samesize))]);
            mymessage = ['Converting Ch ',num2str(currentChannel),' ; File # '];
            mymessage = [mymessage, num2str(fileIter), ' of ', num2str(numFiles),'...'];
            if waitbarFlag == 1, wb = waitbar(0,mymessage); end
            for currentFrame = 1:numFrames;
            fseek(fidr,header_sizes(currentFrame),0);
                %Read the data
                frame_data = fread(fidr,frame_size(1),'int16');
                framePic = reshape(frame_data,xsize,ysize);
                %Write the frame to the 16-bit tiff stack
                if convertToTiffFlag==1,
                    if currentFrame <= numChannels,
                        imwrite(uint16(transpose(framePic)),outputTiffFileName{channel_seq(currentFrame)},...
                            'tiff','Compression','none');
                    else
                        imwrite(uint16(transpose(framePic)),outputTiffFileName{channel_seq(currentFrame)},...
                            'tiff','WriteMode','append','Compression','none');
                    end
                end
                if waitbarFlag == 1, waitbar(currentFrame/numFrames,wb); end
            end %for currentFrame = 1:numFrames,
            if waitbarFlag == 1,close(wb);end
        end %for fileIter = 1:numFiles
end %switch fileType
%close out files and close down figure windows
fclose all;




