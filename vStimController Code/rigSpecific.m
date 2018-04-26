

PSC_moviedirpath = 'c:\VStimConfig\movies\';
PSC_paramdirpath = 'c:\VStimConfig\';
PSC_logdirpath = 'C:\Vstim\VStimLog\';
PSC_DAQ_PC_IP = '132.239.203.6'; %BA
% PSC_moviedirpath = 'C:\SRO DATA\vStim Data\Movies\';
% PSC_paramdirpath = 'C:\SRO DATA\vStim Data\vStimConfig\vStimController\';
% PSC_logdirpath = 'C:\SRO DATA\vStim Data\vStimLog\';
% PSC_DAQ_PC_IP = '132.239.203.99'; %SRO

PSC_REMOTECONTROL_REMOTEPC_IP = '132.239.203.6'; %BA
PSC_REMOTECONTROL_REMOTEPC_PORT = '3458'; %BA
PSC_REMOTECONTROL_LOCALPC_PORT = '3458'; %BA
% PSC_REMOTECONTROL_REMOTEPC_IP = '132.239.203.99'; %SRO
% PSC_REMOTECONTROL_REMOTEPC_PORT = '3458'; %BA
% PSC_REMOTECONTROL_LOCALPC_PORT = '3458'; %BA

% Default screen setup (must specify all 3 fields if VSTIM_RES is defined)
% VSTIM_RES.width = 1680/2;
% VSTIM_RES.height = 1050/2;

% Old KR stim. monitor -- switched on 1/3/2012
% VSTIM_RES.width = 1280;
% VSTIM_RES.height = 800;
% VSTIM_RES.hz = 60;
% VSTIM_RES.screennum = 1;
% New KR stim. monitor -- added on 1/3/2012
% VSTIM_RES.width = 1920;
% VSTIM_RES.height = 1080;

% BEST STIM MONITOR  3/27/13
VSTIM_RES.width = 1440;
VSTIM_RES.height = 900;
VSTIM_RES.hz = 60;
VSTIM_RES.screennum = 2;

% FEF STIM MONITOR 3/27/13
% VSTIM_RES.width=1024;
% VSTIM_RES.height=768;
% VSTIM_RES.hz=60;
% VSTIM_RES.screennum=2;


% VSTIM_RES.width = 1920;
% VSTIM_RES.height = 1200;
%VSTIM_RES.hz = 60;

% Gamma calibration
% if none set to empty. this will use flat look-up table
PSC_GAMMATABLE = '';

% Calibration for vStimController (y = a*x^b, where y luminance cd/m^2 and x
% is pixel value. a and b are fitted constants for a particular monitor)
% a=0.001676; % from Shawn
% b=2.183;
% My calculations from my monitor calibration shown below
% a=1.6462*10^-4; % These are vals for smaller stim. monitor
% b=2.5079;
% a=0.5*(0.00552+0.00684);
% b=0.5*(1.89343+1.84053);

% After calibration on 1/9/2012
a=1.0803;
b=0.9261;

parentfolder(PSC_moviedirpath,1);
parentfolder(PSC_paramdirpath,1);
parentfolder(PSC_logdirpath,1);

% defaults
if ~exist('VSTIM_RES','var'); VSTIM_RES.width = 800; VSTIM_RES.height = 600;  VSTIM_RES.hz = 60; end
