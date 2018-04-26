% Getting camera running

imaqtool
load camera file (Created with Samexpert)
 V = videoinput('coreco',1,'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Dalsa 1m60\dalsa1m60.ccf')
imaqhwinfo(V) % get info

% what acq properties can be controlled?
get(getselectedsource(V))


% getting data
isrunning(V)

triggerconfig(V, 'Manual')
triggerinfo(V) % all possible configs
% more help
http://www.mathworks.com/access/helpdesk/help/toolbox/imaq/index.html?/access/helpdesk/help/toolbox/imaq/bqm3ixt-1.html&http://www.mathworks.com/products/imaq/supportedio13861.html