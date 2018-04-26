bCLEAN60HZ = 0;
for m = 1:length(STF)
    LOADFILE = [PATH STF(m).filename];
  datafile(m).filename = LOADFILE;
    if STF(m).MAXTRIGGER>1;TrigRange = [1 STF(m).MAXTRIGGER];
    else TrigRange = STF(m).MAXTRIGGER; end
   [data dt] = loadDAQData(LOADFILE,LFPchn,TrigRange);
    STF(m).MAXTRIGGER = size(data,2); %update MAXTRIGGER
    
    sAnn = strfind(LOADFILE,'\');sAnn = LOADFILE(sAnn(end)+1:end);
    %%% subsample (Decimate)
    if 1 
        dt = dt*subsample;
        subdata = data(1:subsample:end,:,:);
        clear data
    end
    
    %%% Filter to create LFP
    display(['LOW PASS: ' num2str(dfilter) 'Hz']);
    if bCLEAN60HZ;        dfilter = [dfilter(1),dfilter(2),60 1];         display('Cleaning 60Hz');
    else        dfilter = [dfilter(1),dfilter(2) 0 1];    end
    clear LFP;
    tic
    for i =1 : size(subdata,3);        LFP(:,:,i) = prepdata(squeeze(subdata(:,:,i))',0,dt,dfilter)';    end
    toc
    
    LFPdatafile(m).maxsample{m} = size(LFP,1);
    LFPdatafile(m).LFP= LFP;
    LFPdatafile(m).dfilter = dfilter;
    LFPdatafile(m).filename = datafile(m).filename;
 LFPdatafile(m).LFPchn = LFPchn(m);
 LFPdatafile(m).dt = dt;
end
   
        