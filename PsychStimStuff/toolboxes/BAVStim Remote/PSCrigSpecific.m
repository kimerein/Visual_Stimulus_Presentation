

PSC_moviedirpath = 'c:\VStimConfig\movies\';
PSC_paramdirpath = 'c:\VStimConfig\';
PSC_logdirpath = 'C:\Vstim\VStimLog\';
PSC_DAQ_PC_IP = '132.239.203.6'; %BA

PSC_REMOTECONTROL_REMOTEPC_IP = '132.239.203.6'; %BA
PSC_REMOTECONTROL_REMOTEPC_PORT = '3458'; %BA
PSC_REMOTECONTROL_LOCALPC_PORT = '3458'; %BA

% defualt screen setup (must specify all 3 fields if VSTIM_RES is defined)
% VSTIM_RES.width = 800;
% VSTIM_RES.height = 600;
% VSTIM_RES.hz = 75;
% VSTIM_RES.width = 1920;
% VSTIM_RES.height = 1080;

% BEST STIM MONITOR 3/27/13
VSTIM_RES.width = 1440;
VSTIM_RES.height = 900;
VSTIM_RES.hz = 60;
VSTIM_RES.screennum = 2;

% FEF STIM MONITOR 3/27/13
% VSTIM_RES.width=1024;
% VSTIM_RES.height=768;
% VSTIM_RES.hz=60;
% VSTIM_RES.screennum=2;

% Gamma calibration
% if none set to empty. this will use flat look-up table)
PSC_GAMMATABLE = ''; %E:\Documents and Settings\Bassam\My Documents\Matlab toolboxes\CalibrateGamma\inv_gamma_clut.mat';

parentfolder(PSC_moviedirpath,1);
parentfolder(PSC_paramdirpath,1);
parentfolder(PSC_logdirpath,1);

% defaults
if ~exist('VSTIM_RES','var'); VSTIM_RES.width = 800; VSTIM_RES.height = 600;  VSTIM_RES.hz = 75; end
