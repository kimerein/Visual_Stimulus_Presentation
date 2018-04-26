function  [d si] = axgxload(fn,varargin)
% function  [d si] = axgxload(fn,varargin)
% Desc:
%         loads selected sweeps, channels and time chunk from axgx file
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
 
% Improvements: modify such that import_axgx only loads desired data
% (currently the whole file is loaded and then unwanted data is thrown
% away)


axgxdata = import_axgx(fn);
pvpmod(varargin); % not sure if needed
[d si] = selaxgx(axgxdata,'sweeps',sweeps,'start',start,'stop',stop,'channels',channels);
