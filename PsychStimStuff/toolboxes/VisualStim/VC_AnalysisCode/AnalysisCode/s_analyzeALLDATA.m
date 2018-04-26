% ANALYZED
%%
S_RIGSPECIFIC_SPIKESORT_CONFIG
fid = 1


% spike extraction parameters
Extrparams.invt = 1;
Extrparams.TH_STD = 3.5; % DONT change here!!!
Extrparams.WOI = [0.2 0.6];

% sort data
params.breSort =0;
params.bNoLoad = 0;% don't load previous extracted files


%% 2009-9-03 300 (orientation and spf) not included yet
exptparams{1} = 300;
exptparams{2} = 1;
exptparams{3} = 1;
% Extrparams.invt = -1;
Extrparams.TH_STD = 4.5;

DAQchns = 2;
STF(1).filename = 'Data2009-09-03M1_';
STF(2).filename = 'Data2009-09-03M1_01';
STF(3).filename = 'Data2009-09-03M1_02';
Units = [1]; 
%% 2009-9-16 400 (orientation slightly dirty unit)
exptparams{1} = 400;
exptparams{2} = 1;
exptparams{3} = 1;
% Extrparams.invt = -1;
Extrparams.TH_STD = 5.5;

DAQchns = 2;
STF(1).filename = 'Data2009-09-16M1_05';
STF(2).filename = 'Data2009-09-16M1_06';
STF(3).filename = 'Data2009-09-16M1_07';
Units = [1]; 
%% 2009-9-22 360 not sure about VStim PAram not included yet
exptparams{1} = 360;
exptparams{2} = 1;
exptparams{3} = 1;
% Extrparams.invt = -1;
Extrparams.TH_STD = 5.5;

DAQchns = 2;
STF(1).filename = 'Data2009-09-22M1_02';
Units = [1]; 
%% 2009-9-23 400  (sorted but PROBLEM  VSTim loading) (orientation and spf) not included yet
exptparams{1} = 400;
exptparams{2} = 1;
exptparams{3} = 1;
% Extrparams.invt = -1;
Extrparams.TH_STD = 5.5;

DAQchns = 2;
STF(1).filename = 'Data2009-09-23M1_';
STF(end+1).filename = 'Data2009-09-23M1_03';
STF(end+1).filename = 'Data2009-09-23M1_04';
STF(end+1).filename = 'Data2009-09-23M1_05';
STF(end+1).filename = 'Data2009-09-23M1_06';
Units = [1]; 
%% 2009-9-23 600 (orientation and spf) PART included yet
exptparams{1} = 600;
exptparams{2} = 1;
exptparams{3} = 1;
Extrparams.TH_STD = 5.5;

DAQchns = 2;
STF(1).filename = 'Data2009-09-23M1_08';
STF(end+1).filename = 'Data2009-09-23M1_010';
Units = [1 2]; 
%% 2009-9-23 726 (orientation and spf) not included yet
exptparams{1} = 726;
exptparams{2} = 1;
exptparams{3} = 1;
% Extrparams.invt = -1;
Extrparams.TH_STD = 5.5;

DAQchns = 2;
STF(1).filename = 'Data2009-9-23M1_017';
STF(end+1).filename = 'Data2009-9-23M1_018';
STF(end+1).filename = 'Data2009-9-23M1_019';
Units = [1]; 
%% 2009-10-12 526 (orientation and spf) not included yet
exptparams{1} = 526;
exptparams{2} = 1;
exptparams{3} = 1;
Extrparams.invt = -1;
Extrparams.TH_STD = 5.5;

DAQchns = 2;
STF(1).filename = 'Data2009-10-12M1_';
Units = [1];
%% 2009-10-12 700
exptparams{1} = 700;
exptparams{2} = 1;
exptparams{3} = 1;
% Extrparams.invt = -1;

DAQchns = 2;
STF(1).filename = 'Data2009-10-12M1_03';
Units = [1];
%% 2009-10-22 221 huge number of spikes (very well isolated), no orientation tuning, 
exptparams{1} = 221;
exptparams{2} = 1;
exptparams{3} = 1;
% Extrparams.invt = -1;
Extrparams.TH_STD = 4.5;

DAQchns = 2;
STF(1).filename = 'Data2009-10-22M1_2';
STF(end+1).filename = 'Data2009-10-22M1_3';
STF(end+1).filename = 'Data2009-10-22M1_16';
Units = [1];
%% 2009-10-24 484 (nothing currently included form here ) **  OUTLIERS1, analyzed jsut 7 .. too hard to isolate,
% and low spk rates
exptparams{1} = 484;
exptparams{2} = 1;
exptparams{3} = 1;
% Extrparams.invt = -1;

Extrparams.TH_STD = 5.5;

DAQchns = 2;
STF = [];
% STF(1).filename = 'Data2009-10-24M1_';
% STF(end+1).filename = 'Data2009-10-24M1_05';
STF(end+1).filename = 'Data2009-10-24M1_07';
Units = [1];
%% 2009-10-24 535  (skipped.. come back to) %  WATCH out could be wront STIM cond
exptparams{1} = 535;
exptparams{2} = 1;
exptparams{3} = 1;
% Extrparams.invt = -1;
Extrparams.TH_STD = 5.5;

DAQchns = 2;
STF(1).filename = 'Data2009-10-24M1_08';
STF(end+1).filename = 'Data2009-10-24M1_011';
STF(end+1).filename = 'Data2009-10-24M1_012';
% for i = 31:48 
% STF(end+1).filename = ['Data2009-10-24M1_' num2str(i)];
% end
Units = [1];
cond(1).bmask = 1;
cond(2).bmask = 0;
%% 2009-10-24 550  (Useless not Trigger conditions)
exptparams{1} = 550;
exptparams{2} = 1;
exptparams{3} = 1;
% Extrparams.invt = -1;
Extrparams.TH_STD = 4.5;

DAQchns = 2;
STF(1).filename = 'Data2009-10-24M1_29';
STF(end+1).filename = 'Data2009-10-24M1_30';
for i = 31:46 
STF(end+1).filename = ['Data2009-10-24M1_' num2str(i)];
end
Units = [1];
cond(1).bmask = 1;
cond(2).bmask = 0;
%% 2009-10-24 636  WATCH out could be wront STIM cond (  revisit didn't extract gaussian mask data, and not % sure if 1 or 2 units)
exptparams{1} = 636;
exptparams{2} = 1;
exptparams{3} = 1;
% Extrparams.invt = -1;
Extrparams.TH_STD = 4.5; 
DAQchns = 2;
clear STF; cond = [];
STF(1).filename = 'Data2009-10-24M1_49';
STF(end+1).filename = 'Data2009-10-24M1_50';
STF(end+1).filename = 'Data2009-10-24M1_51';
STF(end+1).filename = 'Data2009-10-24M1_52';
Units = [1 2];
cond(1).bmask = 1;
cond(2).bmask = 0;
%% 2009-10-30 284 
exptparams{1} = 284;
exptparams{2} = 1;
exptparams{3} = 1;
% Extrparams.invt = -1; gd
Extrparams.TH_STD = 4.5;

DAQchns = 2;
clear STF; cond = [];
STF(1).filename = 'Data2009-10-30M1_';
STF(end+1).filename = 'Data2009-10-30M1_02';
STF(end+1).filename = 'Data2009-10-30M1_03';
STF(end+1).filename = 'Data2009-10-30M1_04';
STF(end+1).filename = 'Data2009-10-30M1_07'; % spike may have changed shape , because not
STF(end+1).filename = 'Data2009-10-30M1_05';
STF(end+1).filename = 'Data2009-10-30M1_06';
Units = [1];
cond(1).bmask = 1;
cond(2).bmask = 0;

% MUST SAVE VStim settings ... 
%% 2009-10-30 382 
exptparams{1} = 382;
exptparams{2} = 1;
exptparams{3} = 1;
% Extrparams.invt = -1;
Extrparams.TH_STD = 3; 

DAQchns = 2;
clear STF; cond = [];
STF(1).filename = 'Data2009-10-30M1_C2_01';
STF(end+1).filename = 'Data2009-10-30M1_C2_03';
STF(end+1).filename = 'Data2009-10-30M1_010';
STF(end+1).filename = 'Data2009-10-30M1_13';
STF(end+1).filename = 'Data2009-10-30M1_14';
STF(end+1).filename = 'Data2009-10-30M1_15';
Units = [1];
%% 2009-10-31 337 (non saturating contrast)
exptparams{1} = 337;
exptparams{2} = 1;
exptparams{3} = 1;
% Extrparams.invt = -1;
Extrparams.TH_STD = 3; 

DAQchns = 2;
clear STF; cond = [];
STF(1).filename = 'Data2009-10-31M1_06';
STF(end+1).filename = 'Data2009-10-31M1_07';
STF(end+1).filename = 'Data2009-10-31M1_08';
STF(end+1).filename = 'Data2009-10-31M1_010';
STF(end+1).filename = 'Data2009-10-31M1_011';
Units = [1];
cond(1).bmask = 1;
cond(2).bmask = 0;
%% 2009-10-31 397 (more units but needs more work to extract them)
exptparams{1} = 397;
exptparams{2} = 1;
exptparams{3} = 1;
% Extrparams.invt = -1;

Extrparams.TH_STD = 3; 
clear STF; cond = [];
DAQchns = 2;
STF(1).filename = 'Data2009-10-31M1_012';
Units = [1];
%% 2009-11-02 252  (nothing)
exptparams{1} = 252;
exptparams{2} = 1;
exptparams{3} = 1;
% Extrparams.invt = -1;

DAQchns = 2;
STF(1).filename = 'Data2009-11-02M1_06';
Units = [1];
%% 2009-11-02 225  (nothing)
exptparams{1} = 225;
exptparams{2} = 1;
exptparams{3} = 1;
Extrparams.invt = -1;

DAQchns = 2;
STF(1).filename = 'Data2009-11-02M1_010';
Units = [1];
%% 2009-11-02 162  (nothing) (maybe should be 11 not 011)
exptparams{1} = 162;
exptparams{2} = 1;
exptparams{3} = 1;

DAQchns = 2;
STF(1).filename = 'Data2009-11-02M1_011';
Units = [1];
%% 2009-11-02 276 (nice contrast in varient)
exptparams{1} = 276;
exptparams{2} = 1;
exptparams{3} = 1;
Extrparams.TH_STD = 3; 

clear STF; cond = [];
DAQchns = 2;
STF(1).filename = 'Data2009-11-02M1_15';
STF(2).filename = 'Data2009-11-02M1_16';
STF(3).filename = 'Data2009-11-02M1_17';
STF(4).filename = 'Data2009-11-02M1_18';
Units = [1];
%% 2009-11-02 435 ( not included 1 too few spikes, unit 2 dirty)
exptparams{1} = 435;
exptparams{2} = 1;
exptparams{3} = 1;

Extrparams.TH_STD = 3; 

clear STF; cond = [];
DAQchns = 2;
STF(1).filename = 'Data2009-11-02M1_19';
STF(2).filename = 'Data2009-11-02M1_21';
Units = [1 2];
%% 2009-11-03 203 (orientaiton and spf) (STA)
exptparams{1} = 203;
exptparams{2} = 1;
exptparams{3} = 1;
Extrparams.invt = -1;
Extrparams.TH_STD = 4.5;

clear STF; cond = [];
DAQchns = 2;
STF(1).filename = 'Data2009-11-03M1_';
STF(2).filename = 'Data2009-11-03M1_01';
STF(3).filename = 'Data2009-11-03M1_02';
STF(4).filename = 'Data2009-11-03M1_03';
STF(5).filename = 'Data2009-11-03M1_04';
STF(6).filename = 'Data2009-11-03M1_05';
STF(7).filename = 'Data2009-11-03M1_08';
STF(8).filename = 'Data2009-11-03M1_010';
STF(9).filename = 'Data2009-11-03M1_011';
STF(10).filename = 'Data2009-11-03M1_012';
STF(11).filename = 'Data2009-11-03M1_013';
Units = [1];
%% 2009-11-03 201  (low spike rate, non saturation crf)
exptparams{1} = 201;
exptparams{2} = 1;
exptparams{3} = 1;

Extrparams.invt = 1;
Extrparams.TH_STD = 3.5;

clear STF; cond = [];

DAQchns = 2;
STF(1).filename = 'Data2009-11-03M1_1';
STF(2).filename = 'Data2009-11-03M1_25';
Units = [1];
%% 2009-11-03 228  (
exptparams{1} = 228;
exptparams{2} = 1;
exptparams{3} = 1;

Extrparams.invt = 1;
Extrparams.TH_STD = 3.5;
clear STF; cond = [];

DAQchns = 2;
STF(1).filename = 'Data2009-11-03M1_27';
STF(2).filename = 'Data2009-11-03M1_28';
Units = [2];
%% 2009-11-03 257 
exptparams{1} = 257;
exptparams{2} = 1;
exptparams{3} = 1;

Extrparams.invt = 1;
Extrparams.TH_STD = 3.5;

clear STF; cond = [];

DAQchns = 2;
STF(1).filename = 'Data2009-11-03M1_32';
STF(2).filename = 'Data2009-11-03M1_33';
STF(3).filename = 'Data2009-11-03M1_34';
STF(4).filename = 'Data2009-11-03M1_35';
Units = [1];
%% 2009-11-03 389 
exptparams{1} = 389;
exptparams{2} = 1;
exptparams{3} = 1;

Extrparams.invt = 1;
Extrparams.TH_STD = 3.5;

clear STF; cond = [];

DAQchns = 2;
STF(1).filename = 'Data2009-11-03M1_37';
STF(2).filename = 'Data2009-11-03M1_38';
STF(3).filename = 'Data2009-11-03M1_39';
Units = [1];
%% 2009-11-03 407 
exptparams{1} = 407;
exptparams{2} = 1;
exptparams{3} = 1;

Extrparams.invt = 1;
Extrparams.TH_STD = 3.5;

clear STF; cond = [];

DAQchns = 2;
STF(1).filename = 'Data2009-11-03M1_41';
STF(2).filename = 'Data2009-11-03M1_42';
STF(3).filename = 'Data2009-11-03M1_43';
STF(4).filename = 'Data2009-11-03M1_44';
STF(5).filename = 'Data2009-11-03M1_45';
STF(6).filename = 'Data2009-11-03M1_46';
STF(7).filename = 'Data2009-11-03M1_47';
Units = [1];
%% 2009-11-04 273 
exptparams{1} = 273;
exptparams{2} = 1;
exptparams{3} = 1;
Extrparams.invt = 1;
Extrparams.TH_STD = 3;

clear STF; cond = [];

DAQchns = 2;
STF(1).filename = 'Data2009-11-04M1_06';
STF(2).filename = 'Data2009-11-04M1_07';
STF(3).filename = 'Data2009-11-04M1_08';
STF(4).filename = 'Data2009-11-04M1_09';
Units = [1];
%% 2009-11-06 157 
exptparams{1} = 157;
exptparams{2} = 1;
exptparams{3} = 1;

Extrparams.invt = 1;
Extrparams.TH_STD = 4;
Extrparams.WOI = [0.2 0.4];

clear STF; cond = [];

DAQchns = 2;
STF(1).filename = 'Data2009-11-06M1_02';
STF(2).filename = 'Data2009-11-06M1_03';

Units = [1];
%% 2009-11-06 170 
exptparams{1} = 170;
exptparams{2} = 1;
exptparams{3} = 1;
Extrparams.invt = 1;
Extrparams.TH_STD = 3.5;
Extrparams.WOI = [0.2 0.4];

DAQchns = 2;
STF(1).filename = 'Data2009-11-06M1_05';
STF(2).filename = 'Data2009-11-06M1_06';
STF(3).filename = 'Data2009-11-06M1_07';
Units = [1];
%% 2009-11-09 173 (not included oreintation and spf ) OUTLINER PROBLEM???  with chn 1 
exptparams{1} = 173;
exptparams{2} = 1;
exptparams{3} = 1;

Extrparams.invt = 1;
Extrparams.TH_STD = 4;

Extrparams.WOI = [0.2 0.4];

DAQchns = 2;
% STF(1).filename = 'Data2009-11-09M1_01';
STF(1).filename = 'Data2009-11-09M1_02';
STF(2).filename = 'Data2009-11-09M1_03';
Units = [1 2];
cond = [];
%% 2009-11-09 320 (unit 2 and 3 maybe same unit)
exptparams{1} = 320;
exptparams{2} = 1;
exptparams{3} = 1;

Extrparams.invt = 1;
Extrparams.TH_STD = 4;

Extrparams.WOI = [0.1 0.4];


DAQchns = 2;
STF(1).filename = 'Data2009-11-09M1_05';
STF(2).filename = 'Data2009-11-09M1_06';
Units = [2 3];
%%  2009-11-09 452 % 2 clean units and 1 that contains 2 units together (excuded)
exptparams{1} = 452;
exptparams{2} = 1;
exptparams{3} = 1;


Extrparams.invt = 1;
Extrparams.TH_STD = 4;

Extrparams.WOI = [0.1 0.4];

DAQchns = 2;
STF(1).filename = 'Data2009-11-09M1_08';
STF(2).filename = 'Data2009-11-09M1_09';
STF(3).filename = 'Data2009-11-09M1_010';
STF(4).filename = 'Data2009-11-09M1_011';
STF(5).filename = 'Data2009-11-09M1_012';

% Analysis parameters
Units = [1 2]; 
% cond(1).bmask = 1;
% cond(2).bmask = 0;

% AnalType = 'Orientation'; %'ContrastO'
%%  2009-12-10 230
exptparams{1} = 230;
exptparams{2} = 1;
exptparams{3} = 1;

Extrparams.invt = 1;
Extrparams.TH_STD = 2.5;
clear STF; cond = [];

Extrparams.WOI = [0.2 0.4];


DAQchns = 2;
STF(1).filename = 'Data2009-12-10M1_';

% % Analysis parameters
Units = [1]; 
%%  2009-12-14 352 (broad direction selective unit, and direction sensitive
% unit, .. not clean,(not included) . 
exptparams{1} = 352;
exptparams{2} = 1;
exptparams{3} = 1;

Extrparams.invt = 1;
Extrparams.TH_STD = 4.5;
clear STF; cond = [];

Extrparams.WOI = [0.2 0.4];


DAQchns = 2;
STF(1).filename = 'Data2009-12-14M1_3';
STF(end+1).filename = 'Data2009-12-14M1_4';
STF(end+1).filename = 'Data2009-12-14M1_5';

% % Analysis parameters
Units = [1 4]; 
%%  2009-12-14 406 sta 
exptparams{1} = 406;
exptparams{2} = 1;
exptparams{3} = 1;

Extrparams.invt = 1;
Extrparams.TH_STD = 3.5;
clear STF; cond = [];

Extrparams.WOI = [0.2 0.4];


DAQchns = 2;
STF = [];
for i=6:23
STF(end+1).filename = ['Data2009-12-14M1_' num2str(i)];
end
% % Analysis parameters
Units = [1];
%%  2009-12-14 420
exptparams{1} = 420;
exptparams{2} = 1;
exptparams{3} = 1;

Extrparams.invt = -1;
Extrparams.TH_STD = 3.5;
STF = []; cond = [];
Extrparams.WOI = [0.2 0.4];
DAQchns = 2;
for i=26:28
STF(end+1).filename = ['Data2009-12-14M1_' num2str(i)];
end
% % Analysis parameters
Units = [1];
%%  2009-12-14 511
exptparams{1} = 511;
exptparams{2} = 1;
exptparams{3} = 1;

Extrparams.invt = 1;
Extrparams.TH_STD = 3.5;
 cond = [];
Extrparams.WOI = [0.2 0.4];
DAQchns = 2;
STF = [];
for i=32:33
STF(end+1).filename = ['Data2009-12-14M1_' num2str(i)];
end
% % Analysis parameters
Units = [1:3]
%%  2009-12-14 676 sta
exptparams{1} = 676;
exptparams{2} = 1;
exptparams{3} = 1;

Extrparams.invt = -1;
Extrparams.TH_STD = 3.5;
cond = [];
Extrparams.WOI = [0.2 0.4];
DAQchns = 2;
STF = [];
for i=[40 43:53]
STF(end+1).filename = ['Data2009-12-14M1_' num2str(i)];
end
% % Analysis parameters
Units = [1]; 
%%  2009-12-14 772
exptparams{1} = 772;
exptparams{2} = 1;
exptparams{3} = 1;

Extrparams.invt = -1;
Extrparams.TH_STD = 3.5;
cond = [];
Extrparams.WOI = [0.2 0.4];
DAQchns = 2;
STF = [];
for i=[55 57]
STF(end+1).filename = ['Data2009-12-14M1_' num2str(i)];
end

STF(2).MAXTRIGGER = 76;
% % Analysis parameters
Units = [1]; 
%%
% sexptname = createExptname(STF,exptparams);
FilenameSpkSort = ssortRawData(STF,exptparams,DAQchns,Extrparams,params);
% %  VStim Variables during for each DAQfile recording is saved to a table.
% % VstimTable = updateVStimCondTable(TABLES_SAVEPATH,{STF.filename},0); 
% 
%%  2009-12-20 198 FS cell 
exptparams{1} = 198;
exptparams{2} = 1;
exptparams{3} = 1;

Extrparams.invt = 1;
Extrparams.TH_STD = 4.5;
clear STF; cond = [];

Extrparams.WOI = [0.2 0.4];


DAQchns = 3;
STF(1).filename = 'Data2009-12-20M1_22';
STF(end+1).filename = 'Data2009-12-20M1_23';
STF(end+1).filename = 'Data2009-12-20M1_27';

Units = [1 2 3]
% % Analysis parameters
%%
FilenameSpkSort = ssortRawData(STF,exptparams,DAQchns,Extrparams,params);

%%
analyzeUnitTuningNEW(FilenameSpkSort,STF,Units,'Orientation',cond,0)

%%
analyzeUnitTuningNEW(FilenameSpkSort,STF,Units,'ContrastO',cond,0)
% DO
% analyzeUnitTuningNEW(FilenameSpkSort,STF,Units,'ContrastO',cond,0)
%%
summarizeTuningofUnits(Units,FilenameSpkSort,{STF(2:end).filename},fid)
printOrientfig
%%
xcorrUnits(Units',FilenameSpkSort)
            xcorrmk(PICK_KLUSTERS',PICK_KLUSTERS',stimes,assign,1/dt,'range',200,'bsz',5,'fid',100) ; % 5ms bin
            xcorrmk(PICK_KLUSTERS',PICK_KLUSTERS',stimes,assign,1/dt,'range',20,'fid',10) % 1ms bin
%%
filename = 'UnitTuningStruct';

load(fullfile(TABLES_SAVEPATH,filename));
%%
save(fullfile(TABLES_SAVEPATH,filename),'unitsummary','unitsummaryCond');