%http://www.usd.edu/~schieber/psyc770/USDCalibrate.html
% calibration_phase2.m
%
%CS100A calibration script
%Read-in phase1 luminance data and generate best-fit inverse gamma table
%Save inverse gamma table in phase2_photometry.mat

clear all;
close all;

%load sampled data created via calibration_phase1.m
%indexValues[] and luminanceMeasurements[]
load phase1_photometry.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot sampled luminance values %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf;
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[30 140 800 800]);
subplot(2,2,1);
plot(indexValues,luminanceMeasurements,'+');
hold on;
xlabel('Pixel Values');
ylabel('Luminance (cd/m2)');
strTitle{1}='Sampled Luminance Function';
strTitle{2}='Phase-1 Linear CLUT';
title(strTitle);
axis([0 256 0 max(luminanceMeasurements)]);
axis('square');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate and plot best-fit power function to sampled data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%zero-correct sampled luminance values
lums=luminanceMeasurements-luminanceMeasurements(1);
%normalize sampled luminance values
normalizedLum=lums./max(lums);
%trim zero level
pixels=indexValues(2:end);
normalizedLum=normalizedLum(2:end);

%curve fit empirical luminance values 
fitType = 2;  %extended power function
outputx = [0:255];
[extendedFit,extendedX]=FitGamma(pixels',normalizedLum',outputx',fitType);

%plot sampled luminance and curve fit results
%figure(2);clf;hold on;
subplot(2,2,2); hold on;
plot(pixels,normalizedLum,'+'); %sampled luminance
plot(outputx,extendedFit,'r');  %curve fit results
axis([0 256 0 1]);
xlabel('Pixel Values');
ylabel('Normalized Luminance');
strTitle{1}='Power Function Fit to Sampled Luminance Readings';
strTitle{2}=['Exponent = ',num2str(extendedX(1)),'; Offset = ',num2str(extendedX(2))];
title(strTitle);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate inverse gamma corrected pixel transfer function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pixelMax = max(pixels);
invertedInput=InvertGammaExtP(extendedX,pixelMax,normalizedLum);
%plot inverse gamma function (pixels)
%figure(3); clf; hold on;
subplot(2,2,3); hold on;
plot(pixels,invertedInput,'r+');
axis('square');
axis([0 pixelMax 0 pixelMax]);
plot([0 pixelMax],[0 pixelMax],'r');
xlabel('Pixel Values');
ylabel('Target Pixel Values');
strTitle{1}='Ideal vs. Inverse Gamma Correction';
strTitle{2}=['Exponent = ',num2str(extendedX(1)),'; Offset = ',num2str(extendedX(2))];
title(strTitle);
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate inverse gamma luminance function (based on curve fit above) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxLum = max(luminanceMeasurements);
luminanceRamp=[0:1/255:1];
pow=extendedX(1);
offset=extendedX(2);
invertedRamp=((maxLum-offset)*(luminanceRamp.^(1/pow)))+offset; %invert gamma w/o rounding
%normalize inverse gamma table
invertedRamp=invertedRamp./max(invertedRamp);
%plot inverse gamma function
%figure(4); clf; hold on;
subplot(2,2,4); hold on;
pels=[0:255];
plot(pels,invertedRamp,'r');
axis('square');
axis([0 255 0 1]);
xlabel('Pixel Values');
ylabel('Inverse Gamma Table');
strTitle{1}='Inverse Gamma Table Function';
strTitle{2}=['for Exponent = ',num2str(extendedX(1)),'; Offset = ',num2str(extendedX(2))];
title(strTitle);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% expand inverse gamma to full 3-channel CLUT %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inverseCLUT = repmat(invertedRamp',1,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save inverse gamma table %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save phase2_photometry.mat inverseCLUT

