% Adapted from Izhikevich by Erik Flister
%

clear all
colormap('jet');

% Excitatory (RSi)	Excitatory (RSr)	Inhibitory (FS-horizontal)	Inhibatory (LTS-vertical)
% Excitatory (RSi)	Excitatory (RSr)	Inhibitory (FS-horizontal)	Inhibatory (LTS-vertical)
%-----------------	----------------	--------------------------	-------------------------
Nei=0;   	Ner=0;		Ni=0;			Nlts=2000;		% number of cells

% parameters
% ----------
Cei=100;		Cer=100;		Ci=20;			Clts=100;		% capacitance
Vrei=-60;		Vrer=-60;		Vri=-55;			Vrlts=-50%-60%-50;	% resting voltage
Vtei=-40;		Vter=-40;		Vti=-40;			Vtlts=-42%-40;	% instantaneous threshold
% voltage
kei=.7;		ker=.7;		ki=1;			klts=1;		% k

aei=.03;		aer=.03;					alts=.03;		% a
bei=-2;		ber=5;					blts=8;		% b
cei=-50;		cer=-50;		ci=-45;					% c
dei=100;		der=100;		di=0;			dlts=20;		% d

vpeakei=35;	vpeaker=35;	vpeaki=25;				% v_peak

exiIS=0;		exrIS=0;		inIS=10;			ltsIS=20;		% injected step

numExi=Nei;	numExr=Ner;	numIn=Ni;			numLTS=Nlts;	% num cells to record
% intracellularly
%----------------	---------------	-----------------------------	---------------------------------------------
%    v=v+(I+tau*(k.*(v-Vr).*(v-Vt)-u))./C;

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
    [	-1000	1000	-500	-2000	33
    1000	1000	-500	-2000	33
    150	150	-1200   0	5
    150	150	0	-1000	1]

gapJuncStrengths= ...
    [   0   0   0   0
    0   0   0   0
    0   0   1e-3 0
    0   0   0   0]
%BA test current
testval = 0;
Itest = zeros(1,n);
% Itest(1001:1300) = -200;
% Itest(1:end/4) = testval;
Itest(end/4:1.3*end/4) = testval;
% Itest(end/4:end/2) = testval;
N = total;

%BA
stemp = sprintf('%d %d %d\nltsIS:%d ITest:%d',Nlts,synapseClassStrengths(end-4),synapseClassStrengths(end),ltsIS,testval); % simlabel
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

S=zeros(total,total+1);
G=zeros(total,total);

sideDist=1;
xLocs=sideDist*rand(1,total); % pic xLocation of neuron
yLocs=sideDist*rand(1,total); % pic yLocation of neuron
maxDist=sqrt(2*sideDist^2);
xpLocs=circshift(xLocs,1); % pic xLocation of neuron
ypLocs=circshift(yLocs,1); % pic yLocation of neuron

columnRadius=.05; % in percent of the length of the diagonal of the space
% temp = rand(1,total)*2*pi;
% xTar = 3*columnRadius.*cos(temp);
% yTar =  3*columnRadius.*sin(temp);
% xDists=repmat(xLocs,total,1)-repmat((xLocs+xpLocs)',1,total);
% yDists=repmat(yLocs,total,1)-repmat((yLocs+ypLocs)',1,total);
xDists=repmat(xLocs,total,1)-repmat((xLocs)',1,total);
yDists=repmat(yLocs,total,1)-repmat((yLocs)',1,total);
dists=sqrt(xDists.^2 + yDists.^2)/maxDist; % distance matrix of neurons from all others
% set up gap junctions according to distance
locGWeights = G;
for pre = 1:size(gapJuncStrengths,2)% set up the synapses according to distances
    for post = 1:size(gapJuncStrengths,1)
        locGWeight=zeros(counts(post),counts(pre));
        d=dists((1:counts(post))+sum(counts(1:post-1)),(1:counts(pre))+sum(counts(1:pre-1))); %distances post to pre

        locGWeight = 1-d;
        locGWeight(d>2*(columnRadius)) =0;

        locGWeights(sum(counts(1:post-1))+1:sum(counts(1:post)),...
            sum(counts(1:pre-1))+1:sum(counts(1:pre)))=locGWeight;

        G(sum(counts(1:post-1))+1:sum(counts(1:post)),sum(counts(1:pre-1))+1:sum(counts(1:pre)))= ...
            gapJuncStrengths(post,pre)*rand(counts(post),counts(pre)).*locGWeight;

    end
end

% set up synapses according to distance
locWeights=S;
for pre=1:size(synapseClassStrengths,2)
    for post=1:size(synapseClassStrengths,1)
        if pre==size(synapseClassStrengths,2) %thallumic
            locWeight=ones(counts(post),counts(pre)); %matrix oc local Weights
        else
            locWeight=zeros(counts(post),counts(pre));
            if ~isempty(locWeight)  % Case were some cell types have zero cells
                d=dists((1:counts(post))+sum(counts(1:post-1)),(1:counts(pre))+sum(counts(1:pre-1))); %distances post to pre
                if pre==1
                    %locWeight=(1-d).^6;
                    locWeight=1-d;
                    locWeight(d>(columnRadius))=0;
                elseif pre==2
                    %locWeight=(1-d).^6;
                    locWeight=1-d;
                    locWeight(d>(columnRadius))=0;
                elseif pre==3
                    % 				locWeight=1-d+2*columnRadius;
                    % 				locWeight(d<(2*columnRadius))=0;
                    %                 locWeight(d>(3*columnRadius))=0;
                    %                 locWeight=1-d;
%                     locWeight(d>(2*columnRadius))=0;
                    %                 locWeight = 1;
                    %                
                    locWeight = power(1-d,2);
%                         locWeight(d>(3*columnRadius))=0;
%                       locWeight(d>=(2*columnRadius))=1;
                         locWeight(d<(2*columnRadius))=0;                    locWeight(d>(1*columnRadius))=0;
%                     locWeight(d<=(columnRadius))=1;
%                     locWeight(d<(2*columnRadius))=0;

                    if post ==3
                        temp=1-d;
                        temp(d>columnRadius) = 0;
                        neighbors = find(temp(:,N)>0);
                        temp=1-d;
                        %                    temp(d>(2*columnRadius)) = 0;
                        temp(d<=columnRadius) = 0;
                        temp(d>=2*columnRadius) = 0;

                        extneigh = find(temp(:,N)>0);
                        %                                         neighbors =1;
                        %                                         extneigh = 2;

                    end


                elseif pre==4
                    %                 locWeight(d<(columnRadius))=1;
                    %                     %NEARest NEIGHBOR
                    %                     for jj=1:length(d)
                    %                         temp = min(d(jj,:));
                    %                         ind = find(d(jj,:) == temp(1));
                    %                         locWeight(jj,ind) = 1;
                    %                     end
                    %                         locWeight = 1;
                    %                     temp = max(locWeight);
                    %                     locWeight(locWeight==temp) =1;
                    %                     locWeight(locWeight> temp(1))=0;
                    %% RAdius
                     locWeight = power(1-d,4);
%                     locWeight(end/2:end
                      locWeight(d>=(4*columnRadius))=1;
                                     locWeight(d>(10*columnRadius))=0;
                      locWeight(d<(4*columnRadius))=0;
           
                         locWeight(d<columnRadius)=0;
                       
%                                    locWeight=1; % all to all
                    % BA find cells in within X (2*columnRadius)radius of cell N
                    if post ==4
                        temp=1-d;
                        temp(d>columnRadius) = 0;
                        neighbors = find(temp(:,N)>0);
                        temp=1-d;
                        %                    temp(d>(2*columnRadius)) = 0;
                        temp(d<=columnRadius) = 0;
                        temp(d>=2*columnRadius) = 0;
                        
                        extneigh = find(temp(:,N)>0);
%                                         neighbors =1;
%                                         extneigh = 2;
                        
                    end

                else
                    'error'
                end;
            end

        end;

        locWeights(sum(counts(1:post-1))+1:sum(counts(1:post)),...
            sum(counts(1:pre-1))+1:sum(counts(1:pre)))=locWeight;

        S(sum(counts(1:post-1))+1:sum(counts(1:post)),sum(counts(1:pre-1))+1:sum(counts(1:pre)))= ...
            synapseClassStrengths(post,pre)*rand(counts(post),counts(pre)).*locWeight;
%         S(sum(counts(1:post-1))+1:sum(counts(1:post)),sum(counts(1:pre-1))+1:sum(counts(1:pre)))= ...
%             synapseClassStrengths(post,pre)*1.*locWeight;
        
            %Synaptic weight matrix, randomly connected neurons, locally
        %weighted w/strength of synapseClassStrengths
    end;
end;

newS=S(:,1:end-1);
newS=newS';
posS=max(0,newS);
negS=min(0,newS);
pos=sum(posS); % total ex synapses
neg=sum(negS); % total in synapses

% some broken code that tries to normalize

if false 						% disable synapses
    % so can study response to current steps, etc.
    S=zeros(total,total+1);
    G=zeros(total,total);
end;

close all;
% show example synapses for each cell type
if true
    figure;
    subplot(1,1,1);
    showCell(1+[0 cumsum(counts(1:end-2))],xLocs,yLocs,S);
    pause;
    colormap('jet');
    title(stemp);
end;

if false
    figure;
    subplot(1,3,2);
    imagesc(S);
    title('SynapticWeights')
    subplot(1,3,3);
    imagesc(locWeights);
    title('LocWeights')
    pause;
    % 	subplot(1,1,1);
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
 samps =1; %BA no intra
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

for i=1:n-1	 %each time step

    % decide what's a spike
    ltsPeaks=[];
    if numLTS>0
        ltsPeaks=40-.1*u((Nei+Ner+Ni+1):total);		% lts cells have a funny peak
    end;						% that depends on u – doesn’t affect
    % dynamics, just cosmetic
    %BA THIS is amplitude of spike?
    vpeaks=[vpeakei*ones(Nei,1) ; vpeaker*ones(Ner,1) ; vpeaki*ones(Ni,1) ; ltsPeaks];

    % indices of spikes
    fired=find(v>=vpeaks);

    % pad spikes to peak
    v(fired)=vpeaks(fired);
    Is=sum(S(:,fired),2);				% synaptic current from spiked neurons
    Ig =0;
    if false %Gap junction (not sure it works)
        Ig = sum(G*(repmat(v',total,1) - repmat(v,1,total)))';     %BA?? sum is right?
    end
    if round(i*tau)>currMs
        currMs=round(i*tau);															% only update noise once per ms									% otherwise frequency depends on tau
        % i think we still have a tau dependency
        % for some reason.  ask Izhikevich about this.

        % thalamic input / spontaneous noise
        INoise=S(:,end).*[randn(Nei,1) ; randn(Ner,1) ; randn(Ni,1) ; randn(Nlts,1)];
    end;

    if i>startStim*n && i<stopStim*n
        % BA inject NOISE at
        %         INoise(end-Nlts-1:end) = 0;
        I=INoise+[exiIS*ones(Nei,1) ; exrIS*ones(Ner,1) ;
            inIS*ones(Ni,1) ; ltsIS*ones(Nlts,1)];		% injected input
    else
        I=[zeros(Nei,1) ; zeros(Ner,1) ; zeros(Ni,1) ; zeros(Nlts,1)];
    end;

    I=I+Is+Ig;
    I(fired)=0;					% prevent feedback from autapses
    % ask Izhikevich about how this
    % interacts with tau

    if true %BA inject Itest (ignore all other currents)
        %         I(end-(Nlts+Ni-1)) = I(end-(Nlts+Ni-1)).*0 +Itest(1,i);
        I(neighbors) = I(neighbors).*0 +Itest(1,i);
    end

    % record v, u, and I for posterity
    samp=ceil(i/tausPerRecSamp);
%     if numExi>0
%         exi(:,samp)=max(v(1:numExi),exi(:,samp));
%         exiU(:,samp)=u(1:numExi);
%         exiI(:,samp)=I(1:numExi);
%     end;
%     if numExr>0
%         exr(:,samp)=max(v(Nei+1:Nei+numExr),exr(:,samp));
%         exrU(:,samp)=u(Nei+1:Nei+numExr);
%         exrI(:,samp)=I(Nei+1:Nei+numExr);
%     end;
%     if numIn>0
%         in(:,samp)=max(v((Nei+Ner+1):(Nei+Ner+numIn)),in(:,samp));
%         inU(:,samp)=u((Nei+Ner+1):(Nei+Ner+numIn));
%         inI(:,samp)=I((Nei+Ner+1):(Nei+Ner+numIn));
%         %         inIg(:,samp) = Ig((Nei+Ner+1):(Nei+Ner+numIn));
%     end;
%     if numLTS>0
%         lts(:,samp)=max(v((Nei+Ner+Ni+1):(Nei+Ner+Ni+numLTS)),lts(:,samp));
%         ltsU(:,samp)=u((Nei+Ner+Ni+1):(Nei+Ner+Ni+numLTS));
%         ltsI(:,samp)=I((Nei+Ner+Ni+1):(Nei+Ner+Ni+numLTS));
%         %         ltsIn(:,samp)=INoise((Nei+Ner+Ni+1):(Nei+Ner+Ni+numLTS));
%     end;

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
    v(fired(fired<=Nei+Ner+Ni))=c(fired(fired<=Nei+Ner+Ni)); % treat LTS differently
    v(fired(fired>Nei+Ner+Ni))=-53+.04*u(fired(fired>Nei+Ner+Ni));
%     v(fired(fired>Nei+Ner+Ni))=-40;


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
end;%time steps

% plot output
firings=firings(firings(:,1)>0,:);
size(firings)
% rasterplots
if true
    figure;
    f1=firings(firings(:,2)<=sum(counts(1)),:);
    f2=firings(firings(:,2)<=sum(counts(1:2)) & firings(:,2)>sum(counts(1)),:);
    f3=firings(firings(:,2)<=sum(counts(1:3)) & firings(:,2)>sum(counts(1:2)),:);
    f4=firings(firings(:,2)<=sum(counts(1:4)) & firings(:,2)>sum(counts(1:3)),:);
    %     j=0;j0=0;
    %     for ii=1:size(firings,1)
    %         if (find(firings(ii,2) == neighbors))
    %             j = j +1;   fn(j)=firings(ii,2);
    %         else
    %             j0 = j0+1; fo = firings(ii,2);
    %         end
    %     end
    % by type
%     subplot(1,2,1);
%     hold on
%     plot(f1(:,1),f1(:,2),'b.','MarkerSize',1);
%     plot(f2(:,1),f2(:,2),'k.','MarkerSize',1);
%     plot(f3(:,1),f3(:,2),'g.','MarkerSize',1);
%     plot(f4(:,1),f4(:,2),'r.','MarkerSize',1);
%     hold off
%     xlim([0,T]);
% 
%     % by space
%     subplot(1,2,2);
    hold on
    plot(f4(:,1),xLocs(f4(:,2)),'r.','MarkerSize',10);
    plot(f3(:,1),xLocs(f3(:,2)),'g.','MarkerSize',10);
    plot(f1(:,1),xLocs(f1(:,2)),'b.','MarkerSize',10);
    plot(f2(:,1),xLocs(f2(:,2)),'k.','MarkerSize',10);
    hold off
    xlim([0,T]);

    title(stemp);
end;

% plot u, v, and I
if false
    pause;
    scaleI=10;
    scaleU=25;
    cells=1:1;
    time=1:ceil(((i-1)/tausPerRecSamp));

    % 	subplot(4,1,1);
    % 	plot([exi(cells,time);exiI(cells,time)/scaleI; exiU(cells,time)/scaleU ]');
    %
    % 	subplot(4,1,2);
    % 	plot([exr(cells,time);exrI(cells,time)/scaleI; exrU(cells,time)/scaleU ]');
    %
    % 	subplot(4,1,3);
    plot([in(cells,time);inI(cells,time)/scaleI; inU(cells,time)/scaleU ]');

    % 	subplot(4,1,4);
    % 	plot([lts(cells,time);ltsI(cells,time)/scaleI; ltsU(cells,time)/scaleU ]');
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

    firings(:,1) = firings(:,1)./.15 ;

