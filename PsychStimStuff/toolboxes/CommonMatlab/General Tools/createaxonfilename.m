function [filename] = createaxonfilename(readfileheader,readfilenumber,type)
%function [filename] = createaxonfilename(readfileheader,readfilenumber,[type])
% type = 0 or not included .abf
% type = 1 atf
    if nargin < 3 type = 0; end % default to .abf files
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

    switch type
        case 1
            sExtention= '.atf';
        otherwise
            sExtention = '.abf';
    end
    filename = sprintf('%s%d%s',filename,readfilenumber,sExtention)
