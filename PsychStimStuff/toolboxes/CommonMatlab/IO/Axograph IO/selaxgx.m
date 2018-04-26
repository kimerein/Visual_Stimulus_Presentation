function  [d si] = selaxgx(axgxdata,varargin)
% function  [d si] = selaxgx(axgxdata,varargin)
% Desc:
%         extracts selected sweeps, channels and time chunk from axgx from
%         axgxdata (output of import_axgx.m)
% (modeled input after abfload.m by H. Hentschke)
%                    >>> INPUT VARIABLES >>>
%
% NAME        TYPE, DEFAULT      DESCRIPTION
% fn          char array         abf data file name
% start       scalar, 0          only gap-free-data: start of cutout to be read (unit: sec)
% stop        scalar or char,    only gap-free-data: end of cutout to be read (unit: sec). 
%             'e'                 May be set to 'e' (end of file).
% sweeps      1d-array or char,  only episodic data: sweep numbers to be read. By default, 
%             'a'                 all sweeps will be read ('a')
% channels    cell array         names of channels to be read, like {'IN 0','IN 8'};
%              or char, 'a'       ** make sure spelling is 100% correct (including blanks) **
%                                 if set to 'a', all channels will be read
% OUTPUT: d: timeseries x channel x episode
% sampling rate
si = (axgxdata(1).data(2)-axgxdata(1).data(1))*1e6;% convert from ms to us

scolname = unique({axgxdata.title});
ncol = length(scolname);
% for i = 1:length(scolname)
%     j=1;
%     while j < length(scolname)+1
%         if strfind(axgxdata(j).title,cell2mat(scolname(i)));
%             col(i) = j;
%             j = length(scolname)+1;
%         end
%         j= j+1;
%     end
% end

nsweeps = (length(axgxdata)-1)/(ncol-1);
pvpmod(varargin); % not sure if needed
if exist('sweeps')
    if isstr(sweeps) sweeps = [1:nsweeps];  
    elseif isempty(sweeps); sweeps = nsweeps; end
    nselsweeps = length(sweeps);  
else    nselsweeps = nsweeps;    sweeps = [1:nsweeps]; end

if exist('channels')
    nselcol = length(channels);
    for i = 1:nselcol
        j=1;
        while j < length(scolname)+1
            if strfind(cell2mat(channels(i)),axgxdata(j).title);
                col(i) = j;
                j = length(scolname)+1;
            end
            j= j+1;
        end
    end
else nselcol = ncol -1; col = [2:ncol];end

if ~exist('col')
    stemp = [];
    for i=1:length(axgxdata);        stemp = [stemp ',' axgxdata(i).title];    end
    serr = sprintf('Could not find specified channel name. \n Channel names in file:%s',stemp);
    error(serr)
end
   
if exist('start')% get index   
    if start==0    st_ind = 1;
    elseif isempty(start);  st_ind = 1;
    else st_ind = round(start/(si*1e-6 )); end
else st_ind = 1; end
    
if exist('stop')% get index
    if isstr(stop); sp_ind = size(axgxdata(1).data,1);
    elseif isempty(stop);  sp_ind = size(axgxdata(1).data,1);
    else sp_ind = round(stop/(si*1e-6 )); end
else sp_ind =  size(axgxdata(1).data,1);end

d = zeros([(sp_ind-st_ind)+1 nselcol nselsweeps]);
for i = 1: nselsweeps
    isw = sweeps(i);
    for j = 1:length(col)
        temp = axgxdata(col(j)+(isw-1)*(ncol-1)).data;
        d(:,j,i) = temp(st_ind:sp_ind);
    end
end
