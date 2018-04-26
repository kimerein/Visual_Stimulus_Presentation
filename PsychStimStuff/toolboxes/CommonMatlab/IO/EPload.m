function [d si] = EPload(fn, varargin)
% function EPload(fn, varargin)
% Desc:
%         loads selected sweeps, channels and time chunk from axgx/abf file
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
% BA080507
ext = getfextension(fn);
pvpmod(varargin); % not sure if needed
if ~isempty(strfind(ext,'axg'))||~isempty(strfind(ext,'axd'))
    [d si] = axgxload(fn,'sweeps',sweeps,'start',start,'stop',stop,'channels',channels);
elseif strfind(ext,'abf')
    [d si] = abfload(fn,'sweeps',sweeps,'start',start,'stop',stop,'channels',channels);
else
    error('No data loaded: unknown filename extension');
end