%%%
%%% s_PlotKlusters
% plots like s_plotklusterspikes but used for klusters that have been save
% and are RELOADED (KlusterSpikes_EE5_W7EE5.mat_. (uses same code as s_findFirstSpikeInBurst
triggerEE = 4;

defineDir
writedirheader = DATAANAL_DIR;
writedirheader = [writedirheader sourcefile(1:end-4) '\' 'EE' strrep(num2str(EEgroup),'  ','_') '\'];
if ~isdir(writedirheader)
    mkdir(writedirheader);
end;

breport = 1; bsave = 1;
figIND = 550;
figString = 'Klusters';
nIntracellsignals =1;
clear accept_units; good_type = []; 
for i=1:size(temp_Intra,1)
    accept_units{i} = {i};
end
sunitfilename = 'Kluster';
s_acceptUnits
breport = 1;
s_MSPairFigure
% s_findFirstSpikeInBurst


% Ts =2.0833e+004
% 
% dt =  4.8000e-005