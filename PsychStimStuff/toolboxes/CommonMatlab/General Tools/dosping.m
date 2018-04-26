function status = dosping(sip)
% function status = dosping(sip)
%
% use dos to ping an ip
% sip = string e.g. '132.239.203.4'
%
% status = 1 if packet was recieved
% status = 0 if not
%
% BA 010210
s = sprintf('ping -n 1 %s',sip);
[s r] = dos(s);

if strfind (r,'Lost = 1'); status = 0;
elseif   strfind (r,'Received = 1'); status = 1;
else error('ping may not have worked'); end