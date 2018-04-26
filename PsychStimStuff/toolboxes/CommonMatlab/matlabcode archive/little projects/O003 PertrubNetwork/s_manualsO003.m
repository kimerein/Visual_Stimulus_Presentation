%Enter spiketimes (from pclamp)
spiketimeSw = % [sweepnum timeInsweep ISI]

%
Interval(1)=0
Interval(2) = 0;
[readfileheader readfilenumber]=  readaxonDialog(bSliceCell)
processedFiles = createaxonfilename(readfileheader,readfilenumber(1)) %ATF
N_sweeps = max(spiketimeSw(:,1);
Ts = 1/50e-6;

threshold = 

bpClamp = 1; % set to 1 if data was extracted by pClamp.
%parts of the plot will be left out because data hasn't be loaded into
%Matlab
run sO003_originalPlot
run sO003_PLOT_print