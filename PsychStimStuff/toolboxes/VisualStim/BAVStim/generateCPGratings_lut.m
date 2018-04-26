function [img lut]= generateCPGratings_lut(orient,freq,TempFreq,phase,contrast,duration, degPerPix,sizeX,sizeY, frameRate, black, white,sizeLut)%%% generate images for counterphase gratings, without displaying them%%% based (loosely) on DriftDemo from psychtoolbox%%% cmn 11/05/05% DriftDemo shows a drifting grating.% % See also GratingDemo, NoiseDemo, MovieDemo.% 2/21/02 dgp Wrote it, based on FlickerTest.% 4/23/02 dgp HideCursor & ShowCursor.gray = 0.5*(white+black);if contrast>1    contrast=1;endinc=(white-gray)*contrast;% stimulus parametersframes = duration*frameRate;  % temporal period, in frames, of the drifting gratingframesPerPeriod = frameRate / (TempFreq);wavelength = 1/freq;pixPerDeg = 1/degPerPix;%%% calculate image, a ramp for 0 to 2pi[x,y]=meshgrid(1:sizeX,1:sizeY);angle=orient*pi/180;f= 2*pi/(pixPerDeg*wavelength); % cycles/pixela=cos(angle)*f;b=sin(angle)*f;img = floor(mod(a*x+b*y + phase/(2*pi),2*pi)*sizeLut/(2*pi));%%% calculate lookup table, i.e. cos(kx)sin(wt)ph = (2*pi*(0:(sizeLut-1))/sizeLut)';lut =zeros(sizeLut,3,floor(frames));incph = inc*sin(ph);for i=1:frames    tempphase=(i/framesPerPeriod)*2*pi;    lu = gray + incph*sin(tempphase);  lut(:,:,i) = [lu lu lu];end