function Fdata = convertEdata(Edata)
% function Fdata = convertEdata(Edata)
% Helper function for loading data using loadEPfiles
% takes in Edata (experiment struct) that contains multiple raw data files (.datafile(N)) 
% and converts them to multiple Edata(N) so that loadEPfiles will work
% use with s_Correlate_LFP_pktopk_w_IEI_002.m

for jj = 1:length(Edata.datafile)
    Fdata(jj).sName = Edata.datafile(jj).sName ;
    Fdata(jj).sChan = Edata.datafile(jj).sChan;
    Fdata(jj).sweeps = Edata.datafile(jj).sweeps ;
    % Fdata(jj).sweeps = NaN ;
    Fdata(jj).start = Edata.datafile(jj).start ;
    Fdata(jj).stop = Edata.datafile(jj).stop ;
%     Fdata(jj).swin = Edata.datafile(jj).swin;
end