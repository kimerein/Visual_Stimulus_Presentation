%BA test calibration
% code mostly from http://www.usd.edu/~schieber/psyc770/USDCalibrate.html
% phase2.m
%
% calibration data from useSpyder.m
% spyderCaldata = useSpyder('LCD', 0, [], [])

load('E:\Documents and Settings\Bassam\My Documents\Matlab toolboxes\CalibrateGamma\spyderCaldata.mat');
measurementsPerChannel = 256;
gamma_clut = zeros(256,3);
figure(1);clf;
clear inv_gamma_clut; clear gamma_clut;
for i=1:3
    inds=(1:measurementsPerChannel)+(i-1)*measurementsPerChannel;
    gamma_clut(:,i)= (spyderCaldata(inds,2));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate and plot best-fit power function to sampled data %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %zero-correct sampled luminance values
    lums=gamma_clut(:,i)-gamma_clut(1,i);
    %normalize sampled luminance values
    normalizedLum=lums./max(lums);
    %trim zero level
    pixels=indexValues(2:end);
    normalizedLum=normalizedLum(2:end);
    subplot(1,3,i);plot(normalizedLum,'+');hold on; axis tight;

    fitType = 2;  %extended power function
    outputx = [0:255];
    [extendedFit,extendedX]=FitGamma(pixels',normalizedLum,outputx',fitType);
    plot(outputx,extendedFit,'r'); %curve fit results
    strTitle{1}='Power Function Fit to Sampled Luminance Readings';
    strTitle{2}=['Exponent = ',num2str(extendedX(1)),'; Offset = ',num2str(extendedX(2))];
    title(strTitle);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate inverse gamma luminance function (based on curve fit above) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    maxLum = max(luminanceMeasurements);
    luminanceRamp=[0:1/255:1];
    pow=extendedX(1);
    offset=extendedX(2);
    invertedRamp=((maxLum-offset)*(luminanceRamp.^(1/pow)))+offset; %invert gamma w/o rounding
    %normalize inverse gamma table
    inv_gamma_clut(:,i)=invertedRamp./max(invertedRamp);
    plot(inv_gamma_clut(:,i),'g');
end


%%
% test calibraiton
useSpyder('LCD', 0, inv_gamma_clut, 5)
%%
save('inv_gamma_clut')

