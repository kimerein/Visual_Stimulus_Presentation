dirheader =  'e:\My Documents\Academic\Rotation\Scanziani\Ascii\';
fileheader = '2005_07_12_';
filenumbers = [18];
MAX_CHN = 16;

for (i =1: size(filenumbers,2))
    %Set Filename
    file = strcat(dirheader,fileheader) ;
    if (filenumbers(i)<1000)
        file = strcat(file,'0');
    end
    if (filenumbers(i)<100)
        file = strcat(file,'0');
    end
    if (filenumbers(i)<10)
        file = strcat(file,'0');
    end
    
    %READING
    file = sprintf('%s%d.atf',file,filenumbers(i))
    
    fid = fopen(file,'a');
    temp = textscan(fid,'','headerLines',1,'delimiter',',')
    temp = textscan(fid,'%*f %f','headerLines',1);
    N_col = cell2mat(temp{1,1});
    
    temp = textscan(fid,'"SignalsExported=%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s','headerLines',8,'delimiter',',')
    N_chn = MAX_CHN - sum(cellfun('isempty',temp));
    for i = 1: N_chn
        CHN_name(i) = cell2mat(temp{1,i});
    end    
    N_sweeps = (N_col-1)/N_chn;
    
    temp = textscan(fid,'"SweepStartTimesMS=%s\n','headerLines',7);
    temp = cell2mat(temp{1,1});
    temp = textscan(temp,'%f,');
    Sweep_Times = cell2mat(temp) %row mat
 
    %temp = textscan(fid,'"SweepStartTimesMS=%f\n','headerLines',6,'delimiter',',')
end