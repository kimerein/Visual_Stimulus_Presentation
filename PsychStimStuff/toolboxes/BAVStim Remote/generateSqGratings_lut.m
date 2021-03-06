function [img lut]= generateSqGratings_lut(orient,freq,TempFreq,phase0,contrast,duration, degPerPix,sizeX,sizeY, frameRate, black, white,sizeLut)
%%% generate images for drifting gratings, without displaying them
%%% based (loosely) on DriftDemo from psychtoolbox
%%% cmn 11/05/05

% DriftDemo shows a drifting grating.
% 
% See also GratingDemo, NoiseDemo, MovieDemo.

% 2/21/02 dgp Wrote it, based on FlickerTest.
% 4/23/02 dgp HideCursor & ShowCursor.

% window
gray = 0.5*(white+black);
if contrast>1
    contrast=1;
end
inc=(white-gray)*contrast;

%%% calculate stimulus parameters
frames = duration*frameRate;  % temporal period, in frames, of the drifting grating

%%%framesPerPeriod = frameRate / (TempFreq); use inverse of this to avoid
%%%dividing by zero when tempfreq=0

FrameFreq = TempFreq/frameRate;   %%%grating frequency, in frames;

wavelength = 1/freq;
pixPerDeg = 1/degPerPix;

%%% calculate image, a ramp from 0 to 2pi aligned with grating
[x,y]=meshgrid(1:sizeX,1:sizeY);
angle=orient*pi/180;
angle = pi-angle;  %%% to follow polar coordinate convention


f= 2*pi/(pixPerDeg*wavelength); % cycles/pixel
a=cos(angle)*f;
b=sin(angle)*f;
img = floor(mod(a*x+b*y+phase0*pi/180,2*pi)*sizeLut/(2*pi));

%%% calculate lookup table, i.e. map phase into sinusoid
ph = (2*pi*(0:(sizeLut-1))/sizeLut)';
lut =zeros(sizeLut,3,floor(frames));
for i=1:frames
    phase=(i*FrameFreq)*2*pi;
    temp= ones(size(ph));
   temp(sin(ph+phase)<0)=-1;
 
    lu = gray + inc*temp;
    lut(:,:,i) = [lu lu lu];
end


