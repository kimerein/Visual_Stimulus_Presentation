function savefigure(writedirheader,direc,figdesc,fileformat,fid)
%function savefigure(writedirheader,dirtemp,figdesc,fileformat)
% saves current figure to disk
%     writedirheader: is the directory where the figure should be printed to
%     direc: is the directory where the figure should be printed to (inside
% writedirheader) if it doesn't exist it will be created
%     fileformat: png pdf emf jpg tif
bsvg = 0;


spath = fullfile(writedirheader,direc);
if ~isempty(spath);
    if ~strcmp(spath(end),'\');
        spath = [spath '\'];
    end
end
if (~isdir(spath))
    mkdir(spath);
end
temp = sprintf('%s%s.',spath,figdesc);
ext = fileformat;
temp = [temp ext];

if nargin >= 5 && ~isempty(fid)
    set(0,'CurrentFigure',fid)
else
    fid = gcf;
end

buseprint = 1;
switch(fileformat)
    case 'png'
        s_fileformat = '-dpng';
    case 'pdf'
        s_fileformat = '-dpdf';
        if (length(temp)+3) > 23
            warning('exporting to pdf does not seem to take files longer then 23 characters')
        end
    case 'emf'
                s_fileformat = '-dmeta';
    case 'jpg'
                s_fileformat = '-djpeg';
    case 'tif'
                s_fileformat = '-dtiffn';
    case 'eps'
%                 s_fileformat = '-depsc2';
        s_fileformat = 'eps';
        buseprint=0;
        exportfig(fid,temp,'Format',s_fileformat,'color','cmyk','Height',8.5,'Width',11);

    case 'svg'
        plot2svg(temp);
        buseprint = 0;
        
    case 'fig'
        buseprint=0;
        P = get(gcf,'Position');
        set(gcf,'Position',get( 0, 'ScreenSize' ));
        saveas(fid,temp);
        set(gcf,'Position',P);

end

if buseprint
    print(s_fileformat,temp);
end

