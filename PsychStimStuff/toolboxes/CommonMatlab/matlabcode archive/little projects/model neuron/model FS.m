Ni = 1 ; index_FS = 1;%index of FS cells in matrix of all cells
%-----------------	----------------	--------------------------	-------------------------
Nei=1;   	Ner=0;		Ni=1;			Nlts=0;		% number of cells

% parameters
% ----------
Cei=100;		Cer=100;		Ci=20;             Clts=100;		% capacitance
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
%--------------------
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

	vpeaks=[vpeakei*ones(Nei,1) ; vpeaker*ones(Ner,1) ; vpeaki*ones(Ni,1) ; ltsPeaks];

% preallocate intracellular records

T=500; tau=.1;			% time span and step (ms)
n = round(T/tau);

stimprot = repmat([500:-100:-400]',1,n);
nsw = size(stimprot,1);
IStim= repmat([zeros(1,n/5) ones(1,n*3/5) zeros(1,n/5)],nsw,1).*stimprot;
%%
vsave = zeros(n,nsw,2);
for j=1:nsw % for several stimuli
    for i=1:n	 %each time step
        vsave(i,j) = v;
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
        v(fired)=vpeaks(fired);

        INoise = 0;randn(1,1)
        I=INoise+IStim(j,n);

        u(fired)=u(fired)+d(fired);

        v(fired(fired<=Nei+Ner+Ni))=c(fired(fired<=Nei+Ner+Ni)); % treat LTS differently
        v(fired(fired>Nei+Ner+Ni))=-53+.04*u(fired(fired>Nei+Ner+Ni));

        % Euler method
        v=v+(I+tau*(k.*(v-Vr).*(v-Vt)-u))./C;

        normals=[1:(Nei+Ner) (Nei+Ner+Ni+1):total];
        u(normals)=u(normals)+tau*a(normals).*(b(normals).*(v(normals)-Vr(normals)) - u(normals));

        % deal with FS cells seperately
        w = v(index_FS);
        U=zeros(Ni,1);
        Vb=-55;
        V=find(w>=Vb);
        U(V)=.025*(w(V)-Vb).^3;
        if Ni>0
            u(index_FS)=u(index_FS)+tau*.2*(U-u(index_FS));
        end;
        if any(any(isnan([u v]))) || any(any(isinf([u v])))
            'exploded!'
            break;
        end;

    end
end