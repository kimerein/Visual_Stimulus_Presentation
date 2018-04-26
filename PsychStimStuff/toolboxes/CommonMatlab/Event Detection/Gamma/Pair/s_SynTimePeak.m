%% PAirs during gamma
%% 
%  TIME OF PEAK ex - IN

breport =1;
defineDir
expttype = 'KinateOsc';
writedirpath = [DATAANAL_DIR expttype '\'];
writedirheader = [DATAANAL_DIR expttype '\'];


A =1;
B = 2;
w = 4;
%% make event time monotonicsly increasing rather then restart at each
%% sweep
lengthsw = 4000 ;%% USER ENTRY %% this should be the length of a sweep in ms
temp1 = event(A).time+(event(A).Trace-1).*lengthsw;
temp2 = event(B).time+(event(B).Trace-1).*lengthsw;
ind = findSimEvent(temp1,temp2,w);
N= length(ind);
bp         fslot=1;

if event(A).V < -40 && event(B).V > -40
    deltatime = event(B).time(ind(:,2)) -event(A).time(ind(:,1))
elseif  event(B).V < -40 && event(A).V > -40
    deltatime = event(A).time(ind(:,1)) -event(B).time(ind(:,2))
else
    display('Not In and Ex')
    bplot = 0;
end
if bplot
    figure(20)
    hist(deltatime,40);
    title(['m= ' num2str(mean(deltatime),'%1.1f') ' s= ' num2str(std(deltatime),'%1.1f') ])
    savetype = 'emf';
    figString = ['Peaktime_' event(A).expnum '_' event(A).sfilename '_' event(B).sfilename ]
    if breport
        dirtemp = 'Pairs';
        figdesc = [figString];
        savefigure(writedirheader,dirtemp,figdesc,savetype)
    end
end
