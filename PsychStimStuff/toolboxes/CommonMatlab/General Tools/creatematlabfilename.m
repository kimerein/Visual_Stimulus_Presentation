function [filename] = creatematlabfilename(readfileheader,readfilenumber)
%Set Filename
        if ~isempty(regexp(readfileheader,'_'))
           if (readfilenumber<1000)
            filename = strcat(readfileheader,'0');
           end
        
        else
          filename = readfileheader;

        end
        if (readfilenumber<100)
            filename = strcat(filename,'0');
        end
        if (readfilenumber<10)
            filename = strcat(filename,'0');
        end
        filename = sprintf('%s%d.mat',filename,readfilenumber);
