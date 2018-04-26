% Adapted from Izhikevich by Erik Flister
%

clear all
colormap('jet');

% Excitatory (RSi)	Excitatory (RSr)	Inhibitory (FS-horizontal)	Inhibatory (LTS-vertical)
%-----------------	----------------	--------------------------	-------------------------
Nei=0;   	Ner=0;		Ni=0;			Nlts=500;		% number of cells

% parameters
% ----------
Cei=100;		Cer=100;		Ci=20;			Clts=100;		% capacitance
Vrei=-60;		Vrer=-60;		Vri=-55;			Vrlts=-56;	% resting voltage
Vtei=-40;		Vter=-40;		Vti=-40;			Vtlts=-42;	% instantaneous threshold
% voltage
kei=.7;		ker=.7;		ki=1;			klts=1;		% k

aei=.03;		aer=.03;					alts=.03;		% a
bei=-2;		ber=5;					blts=8;		% b
cei=-50;		cer=-50;		ci=-45;					% c
dei=100;		der=100;		di=0;			dlts=20;		% d

vpeakei=35;	vpeaker=35;	vpeaki=25;				% v_peak

exiIS=0;		exrIS=0;		inIS=0;			ltsIS=0;		% injected step

numExi=Nei;	numExr=Ner;	numIn=Ni;			numLTS=Nlts;	% num cells to record
% intracellularly
%----------------	---------------	-----------------------------	---------------------------------------------

T=500; tau=.15;			% time span and step (ms)
tausPerRecSamp=1;			% how often to record data (this many steps per sample)
maxExpectedSpikeRate=100;		% in Hz (limits number of spikes we can record)
startStim=.01; stopStim=.85;		% percent range through T for which noise and step stim are active

counts=[Nei, Ner, Ni, Nlts, 1];
typeinds=cumsum(counts);
total=Nei+Ner+Ni+Nlts;
expectedNumSpikes=total*maxExpectedSpikeRate*T/1000;
firings=zeros(expectedNumSpikes,2); 					% preallocate extracellular records
n=round(T/tau);                         				% number of simulation steps

synapseClassStrengths= ...
    [	0	0	0	-0	0
    0	0	0	-0	0
    0	0	0	0	0
    0	0	0	-0	50]

% synapses:
% every time a neuron fires, it gives a instantaneous change to the v of each of the other neurons
% according to the entries in its column in S.  these matricies have the following structure:
%
%	pre: 	ei's	er's	ih's	iv's	thalamus (external excitatory drive or spontaneous noise)
% post:
% ei's
% er's
% ih's
% iv's

%BA test current
Itest = zeros(1,n);
% Itest(1001:1300) = -200;
Itest(101:1000) = 40;

S=zeros(total,total+1);

sideDist=1;
xLocs=sideDist*rand(1,total);
yLocs=sideDist*rand(1,total);
maxDist=sqrt(2*sideDist^2);

xDists=repmat(xLocs,total,1)-repmat(xLocs',1,total);
yDists=repmat(yLocs,total,1)-repmat(yLocs',1,total);
dists=sqrt(xDists.^2 + yDists.^2)/maxDist;

columnRadius=.05; % in percent of the length of the diagonal of the space

% set up the synapses according to distances
locWeights=S;
for pre=1:size(synapseClassStrengths,2)
    for post=1:size(synapseClassStrengths,1)
        if pre==size(synapseClassStrengths,2)
            locWeight=ones(counts(post),counts(pre));
        else
            locWeight=zeros(counts(post),counts(pre));
            d=dists((1:counts(post))+sum(counts(1:post-1)),(1:counts(pre))+sum(counts(1:pre-1)));
            if pre==1
                %locWeight=(1-d).^6;
                locWeight=1-d;
                locWeight(d>(columnRadius))=0;
            elseif pre==2
                %locWeight=(1-d).^6;
                locWeight=1-d;
                locWeight(d>(columnRadius))=0;
            elseif pre==3
                locWeight=1-d+2*columnRadius;
                locWeight(d<(2*columnRadius))=0;
                locWeight(d>(3*columnRadius))=0;
                locWeight=1;
            elseif pre==4
                locWeight(d<(columnRadius))=1;
                 locWeight=1;
           else
                'error'
            end;

        end;

        locWeights(sum(counts(1:post-1))+1:sum(counts(1:post)),...
            sum(counts(1:pre-1))+1:sum(counts(1:pre)))=locWeight;

        S(sum(counts(1:post-1))+1:sum(counts(1:post)),sum(counts(1:pre-1))+1:sum(counts(1:pre)))= ...
            synapseClassStrengths(post,pre)*rand(counts(post),counts(pre)).*locWeight;
    end;
end;

newS=S(:,1:end-1);
newS=newS';
posS=max(0,newS);
negS=min(0,newS);
pos=sum(posS);
neg=sum(negS);

% some broken code that tries to normalize
% synaptic weights to reasonable values
if false
    avgMax=75;
    maxpos=avgMax*sum(counts(1:2));
    inToEx=.3;

    posS=posS/(max(pos)/maxpos);
    pos=sum(posS);
    avgpos=mean(pos);
    avgneg=-1*inToEx*mean(pos);
    negS=negS/(mean(neg)/avgneg);
    neg=sum(negS);

    newS=(posS'+negS');
    S(:,1:end-1)=newS;
end;

if false 						% disable synapses
    % so can study response to current steps, etc.
    S=zeros(total,total+1);
end;

% show histograms of synaptic weights
if false
    subplot(3,1,1);
    hist([pos(1:typeinds(1)); pos(typeinds(1)+1:typeinds(2));
        [NaN*ones(1,counts(1)-counts(3)), pos(typeinds(2)+1:typeinds(3))];
        [NaN*ones(1,counts(1)-counts(4)), pos(typeinds(3)+1:typeinds(4))]]',100);
    subplot(3,1,2);
    hist([neg(1:typeinds(1)); neg(typeinds(1)+1:typeinds(2));
        [NaN*ones(1,counts(1)-counts(3)), neg(typeinds(2)+1:typeinds(3))];
        [NaN*ones(1,counts(1)-counts(4)), neg(typeinds(3)+1:typeinds(4))]]',100);
    subplot(3,1,3);
    hist([reshape(newS(:,1:typeinds(1)),1,total*counts(1)) ;
        reshape(newS(:,typeinds(1)+1:typeinds(2)),1,total*counts(2)) ;
        [reshape(newS(:,typeinds(2)+1:typeinds(3)),1,total*counts(3)) ,...
        NaN*ones(1,total*(counts(1)-counts(3)))];
        [reshape(newS(:,typeinds(3)+1:typeinds(4)),1,total*counts(4)) ,...
        NaN*ones(1,total*(counts(1)-counts(4)))]]',1000);
    pause;
end;

% show example synapses
if false
    subplot(1,1,1);
    showCell(1+[0 cumsum(counts(1:end-2))],xLocs,yLocs,S);
    pause;
    colormap('jet');
end;

if false
    subplot(1,3,2);
    imagesc(S);
    subplot(1,3,3);
    imagesc(locWeights);
    pause;
    subplot(1,1,1);
end;

v=[Vrei*ones(Nei,1); Vrer*ones(Ner,1); Vri*ones(Ni,1); Vrlts*ones(Nlts,1)];			% initial v
u=[zeros(Nei,1); zeros(Ner,1); zeros(Ni,1); zeros(Nlts,1)];					% initial u

C=[Cei*ones(Nei,1); Cer*ones(Ner,1); Ci*ones(Ni,1); Clts*ones(Nlts,1)];
k=[kei*ones(Nei,1); ker*ones(Ner,1); ki*ones(Ni,1); klts*ones(Nlts,1)];
Vr=[Vrei*ones(Nei,1); Vrer*ones(Ner,1); Vri*ones(Ni,1); Vrlts*ones(Nlts,1)];
Vt=[Vtei*ones(Nei,1); Vter*ones(Ner,1); Vti*ones(Ni,1); Vtlts*ones(Nlts,1)];

a=[aei*ones(Nei,1); aer*ones(Ner,1); zeros(Ni,1); alts*ones(Nlts,1)];
b=[bei*ones(Nei,1); ber*ones(Ner,1); zeros(Ni,1); blts*ones(Nlts,1)];
c=[cei*ones(Nei,1); cer*ones(Ner,1); ci*ones(Ni,1); zeros(Nlts,1)];
d=[dei*ones(Nei,1); der*ones(Ner,1); di*ones(Ni,1); dlts*ones(Nlts,1)];

% preallocate intracellular records
largeNeg=-99999999;
samps=ceil(n/tausPerRecSamp);
exi=largeNeg*ones(numExi,samps);
exiU=exi;
exiI=exi;
exr=largeNeg*ones(numExr,samps);
exrU=exr;
exrI=exr;
in=largeNeg*ones(numIn,samps);
inU=in;
inI=in;
lts=largeNeg*ones(numLTS,samps);
ltsU=lts;
ltsI=lts;

ltsResetViolations=0;

spikeNum=1;
currMs=-1;
offsetTimes=[0];

for i=1:n-1

    % decide what's a spike
    ltsPeaks=[];
    if numLTS>0
        ltsPeaks=40-.1*u((Nei+Ner+Ni+1):total);		% lts cells have a funny peak
    end;						% that depends on u – doesn’t affect
    % dynamics, just cosmetic

    vpeaks=[vpeakei*ones(Nei,1) ; vpeaker*ones(Ner,1) ; vpeaki*ones(Ni,1) ; ltsPeaks];

    % indices of spikes
    fired=find(v>=vpeaks);

    % pad spikes to peak
    v(fired)=vpeaks(fired);
    Is=sum(S(:,fired),2);				% synaptic current from spiked neurons

    if round(i*tau)>currMs
        currMs=round(i*tau);															% only update noise once per ms									% otherwise frequency depends on tau
        % i think we still have a tau dependency
        % for some reason.  ask Izhikevich about this.

        % thalamic input / spontaneous noise
        INoise=S(:,end).*[randn(Nei,1) ; randn(Ner,1) ; randn(Ni,1) ; randn(Nlts,1)];

        if false					% noise test
            numOffsets=7;
            numScales=11;

            maxScale=5;
            minScale=0;

            maxOffset=30;
            minOffset=-30;

            offsets=minOffset:((maxOffset-minOffset)/(numOffsets-1)):maxOffset;
            scales=minScale:((maxScale-minScale)/(numScales-1)):maxScale;

            available=stopStim*n-startStim*n;
            segLength=floor(available/(numOffsets*numScales));
            segTimes=ceil(startStim*n):segLength:floor(stopStim*n);
            currSeg=max(find(segTimes<i));
            if ~isempty(currSeg)
                if currSeg>=1 && currSeg<=numOffsets*numScales
                    offsetNum=ceil(currSeg/(numScales));
                    offset=offsets(offsetNum);
                    noiseScale=scales(mod(currSeg-1,numScales)+1);
                    INoise=offset+noiseScale*S(:,end).*[randn(Nei,1) ;
                        randn(Ner,1) ;
                        randn(Ni,1) ;
                        randn(Nlts,1)];
                    if currSeg==1 && offsetTimes(1)==0;
                        offsetTimes(1)=i;
                    end;
                    offsetTimes(offsetNum+1)=i;
                else
                    INoise=0;
                end;
            end;
        end;
    end;

    if i>startStim*n && i<stopStim*n
        I=INoise+[exiIS*ones(Nei,1) ; exrIS*ones(Ner,1) ;
            inIS*ones(Ni,1) ; ltsIS*ones(Nlts,1)];		% injected input
    else
        I=[zeros(Nei,1) ; zeros(Ner,1) ; zeros(Ni,1) ; zeros(Nlts,1)];
    end;

    I=I+Is;
    I(fired)=0;					% prevent feedback from autapses
    % ask Izhikevich about how this
    % interacts with tau
    if true %BA
        I = I.*0 + Itest(i);
    end

    % record v, u, and I for posterity
    samp=ceil(i/tausPerRecSamp);
    % 	if numExi>0
    % 		exi(:,samp)=max(v(1:numExi),exi(:,samp));
    % 		exiU(:,samp)=u(1:numExi);
    % 		exiI(:,samp)=I(1:numExi);
    % 	end;
    % 	if numExr>0
    % 		exr(:,samp)=max(v(Nei+1:Nei+numExr),exr(:,samp));
    % 		exrU(:,samp)=u(Nei+1:Nei+numExr);
    % 		exrI(:,samp)=I(Nei+1:Nei+numExr);
    % 	end;
    % 	if numIn>0
    % 		in(:,samp)=max(v((Nei+Ner+1):(Nei+Ner+numIn)),in(:,samp));
    % 		inU(:,samp)=u((Nei+Ner+1):(Nei+Ner+numIn));
    % 		inI(:,samp)=I((Nei+Ner+1):(Nei+Ner+numIn));
    % 	end;
    % 	if numLTS>0
    % 		lts(:,samp)=max(v((Nei+Ner+Ni+1):(Nei+Ner+Ni+numLTS)),lts(:,samp));
    % 		ltsU(:,samp)=u((Nei+Ner+Ni+1):(Nei+Ner+Ni+numLTS));
    % 		ltsI(:,samp)=I((Nei+Ner+Ni+1):(Nei+Ner+Ni+numLTS));
    % 	end;

    % indicate progress
    pct=100*i/(T/tau);
    if rand>.98
        disp(sprintf('%3.1f\t%d\t%d\t%3.1f\t%d',...
            pct, currMs, spikeNum-1,(spikeNum-1)/(currMs/1000)/total ,ltsResetViolations));
    end;

    % check that we haven't exceeded
    % spike allocation and record spike times
    if spikeNum+length(fired)<=expectedNumSpikes
        firings(spikeNum:(spikeNum+length(fired)-1),:)=[i*tau+0*fired,fired];
        spikeNum=spikeNum+length(fired);
    else
        'too many spikes -- aborting simulation at'
        i
        spikeNum+length(fired)
        expectedNumSpikes
        break;
    end

    % update v and u for spiked neurons
    v(fired(fired<=Nei+Ner+Ni))=c(fired(fired<=Nei+Ner+Ni));
    v(fired(fired>Nei+Ner+Ni))=-53+.04*u(fired(fired>Nei+Ner+Ni));

    limit=-40;
    if any(v(fired)>=limit)				% the lts cells have a funny v reset that
        % depends on u.  sometimes it gets too high,
        % even above the spike detecting threshold
        % so we need to limit it to a reasonable value.
        % ask Izhikevich about fixing this.
        cell=find(v(fired)>=limit);
        cell=fired(cell);
        if any(cell<=sum(counts(1:3)))
            'it wasnt an lts cell'
            cell
            break;
        end;
        v(cell)=limit;
        ltsResetViolations=ltsResetViolations+size(cell,1);
    end;

    u(fired)=u(fired)+d(fired);

    % forward Euler method
    v=v+(I+tau*(k.*(v-Vr).*(v-Vt)-u))./C;

    normals=[1:(Nei+Ner) (Nei+Ner+Ni+1):total];
    u(normals)=u(normals)+tau*a(normals).*(b(normals).*(v(normals)-Vr(normals)) - u(normals));

    weirds=(Nei+Ner+1):(Nei+Ner+Ni);			% FS cells need this alteration of the model
    w=v(weirds);
    U=zeros(Ni,1);
    Vb=-55;
    V=find(w>=Vb);
    U(V)=.025*(w(V)-Vb).^3;
    if Ni>0
        u(weirds)=u(weirds)+tau*.2*(U-u(weirds));
    end;

    % check for explosion
    if any(any(isnan([u v]))) || any(any(isinf([u v])))
        'exploded!'
        break;
    end;

    %fflush(stdout);
end;

% plot output
firings=firings(firings(:,1)>0,:);
size(firings)

% rasterplots
if false
    f1=firings(firings(:,2)<=sum(counts(1)),:);
    f2=firings(firings(:,2)<=sum(counts(1:2)) & firings(:,2)>sum(counts(1)),:);
    f3=firings(firings(:,2)<=sum(counts(1:3)) & firings(:,2)>sum(counts(1:2)),:);
    f4=firings(firings(:,2)<=sum(counts(1:4)) & firings(:,2)>sum(counts(1:3)),:);

    % by type
    subplot(1,2,1);
    hold on
    plot(f1(:,1),f1(:,2),'b.','MarkerSize',1);
    plot(f2(:,1),f2(:,2),'k.','MarkerSize',1);
    plot(f3(:,1),f3(:,2),'g.','MarkerSize',1);
    plot(f4(:,1),f4(:,2),'r.','MarkerSize',1);
    hold off

    % by space
    subplot(1,2,2);
    hold on
    plot(f4(:,1),xLocs(f4(:,2)),'r.','MarkerSize',1);
    plot(f3(:,1),xLocs(f3(:,2)),'g.','MarkerSize',1);
    plot(f1(:,1),xLocs(f1(:,2)),'b.','MarkerSize',1);
    plot(f2(:,1),xLocs(f2(:,2)),'k.','MarkerSize',1);
    hold off
end;

% plot u, v, and I
if false
    pause;
    scaleI=10;
    scaleU=25;
    cells=1:1;
    time=1:ceil(((i-1)/tausPerRecSamp));

    subplot(4,1,1);
    plot([exi(cells,time);exiI(cells,time)/scaleI; exiU(cells,time)/scaleU ]');

    subplot(4,1,2);
    plot([exr(cells,time);exrI(cells,time)/scaleI; exrU(cells,time)/scaleU ]');

    subplot(4,1,3);
    plot([in(cells,time);inI(cells,time)/scaleI; inU(cells,time)/scaleU ]');

    subplot(4,1,4);
    plot([lts(cells,time);ltsI(cells,time)/scaleI; ltsU(cells,time)/scaleU ]');
end;

% traces of i and v and histogram of Vm
if false
    pause;
    cells=1:1;
    scaleI=2;
    maxtime=ceil(((i-1)/tausPerRecSamp));
    time=1:maxtime;

    subplot(4,2,1);
    plot([exiI(cells,time)/scaleI;exi(cells,time)]');

    subplot(4,2,2);
    hist(reshape(exi(:,time),Nei*maxtime,1),100);

    subplot(4,2,3);
    plot([exrI(cells,time)/scaleI;exr(cells,time)]');

    subplot(4,2,4);
    hist(reshape(exr(:,time),Ner*maxtime,1),100);

    subplot(4,2,5);
    plot([inI(cells,time)/scaleI;in(cells,time)]');

    subplot(4,2,6);
    hist(reshape(in(:,time),Ni*maxtime,1),100);

    subplot(4,2,7);
    plot([ltsI(cells,time)/scaleI;lts(cells,time)]');

    subplot(4,2,8);
    hist(reshape(lts(:,time),Nlts*maxtime,1),100);
end;

% plot results of noise test
% you have to set A and AI to the cell you want (i.e., exi and exiI)
if false
    A=exi;
    AI=exiI;

    cells=1:1;
    scaleI=5;
    maxtime=ceil(((i-1)/tausPerRecSamp));
    time=1:maxtime;

    subplot(2,1,1);
    plot([AI(cells,time)/scaleI;A(cells,time);zeros(1,maxtime)]');

    subplot(2,1,2);
    thisA=NaN*zeros(size(A,1)*max(diff(offsetTimes)),length(offsetTimes)-1);
    for st=1:length(offsetTimes)-1
        rows=Nei*(offsetTimes(st+1)-offsetTimes(st));
        thisA(1:rows,st)=reshape(A(:,offsetTimes(st):offsetTimes(st+1)-1),rows,1);
    end;
    hist(thisA,100);
end;

%extracellular movie
if false
    subplot(1,1,1);
    maxtime=ceil(((i-1)/tausPerRecSamp));
    scale=300;
    frame=1;

    done=0;
    for j=1:5:maxtime
        hold off
        plot(scale*xLocs,scale*yLocs,'k.','MarkerSize',10);
        hold on
        plot(scale*xLocs(1:sum(counts(1:2))),scale*yLocs(1:sum(counts(1:2))),'b.','MarkerSize',10);
        fires=firings(:,1);
        fires=fires<j&fires>done;
        plot(scale*xLocs(firings(fires,2)),scale*yLocs(firings(fires,2)),'g.','MarkerSize',10);
        fires=firings(fires,2);
        fires=fires(fires>sum(counts(1:2)));
        plot(scale*xLocs(fires),scale*yLocs(fires),'r.','MarkerSize',10);
        done=j;
        F(frame)=getframe;
        frame=frame+1;
        j
    end;
    movie(F);
    movie2avi(F,'extracellular.avi');
end;

%intracellular movie
if false
    subplot(1,1,1);
    %cells are rows, time is columns
    samps=ceil(((i-1)/tausPerRecSamp));
    recording=[exi(:,1:samps);exr(:,1:samps);in(:,1:samps);lts(:,1:samps)];
    scale=300;
    frame=1;

    lowest=min(min(recording));
    highest=max(max(recording));
    for k=1000:5:size(recording,2)
        a=lowest*ones(scale);

        for j=1:size(recording,1)
            a(ceil(scale*xLocs(j)),ceil(scale*yLocs(j)))=recording(j,k);
        end;
        a(1,1)=highest;
        imagesc(a);
        F(frame)=getframe;
        frame=frame+1;
        k
    end;
    movie(F);
    movie2avi(F,'intracellular.avi');
end;



