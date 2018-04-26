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
   % fclose(fid);
    file = sprintf('%s%d.atf',file,filenumbers(i))
    
    fid = fopen(file,'r');
    temp = textscan(fid,'%*f %f\n','headerLines',1);
    N_col = cell2mat(temp);
    temp = textscan(fid,'"SignalsExported=','headerLines',7)
    temp = textscan(fid,'%s""') %if I don't make it into a string the cell have '' rather then are empty... (dont' know why)
    temp = textscan(cell2mat(temp{1,1}),'%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s','delimiter',',')
    N_chn = MAX_CHN - sum(cellfun('isempty',temp));
    Chn_name = char(zeros(N_chn,6));
    for ii = 1: N_chn
        a = cell2mat(temp{1,ii})
        Chn_name(ii,1:6) = a(1,1:6)  %how do I copy strings
    end    
    N_sweeps = (N_col-1)/N_chn
  %need to go back to the beginning of the file
    fid = fopen(file,'r');
    temp = textscan(fid,'"SweepStartTimesMS=%s\n','headerLines',7)
    str_t = cell2mat(temp{1,1});
    temp = textscan(str_t,'%f,');
    Sweep_Times = cell2mat(temp); %row mat
 
end