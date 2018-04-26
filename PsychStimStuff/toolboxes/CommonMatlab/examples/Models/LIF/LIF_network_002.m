% inhibitory and excitatory neurons
% with synaptic Conductance with time constant
% BA092308

Ne = 1;         Ni = 1;

Vte=-45e-3;     Vti=-45e-3;	% Threshold
Vre=-44e-3;     Vri=-50e-3; % Resting Potential
ce = -55e-3;    ci = -55e-3; % reset potential
Re = 100e6;     Ri = 20e6; % Ohms
Ce = 0.1e-9;    Ci = 0.1e-9; %Farads
Vie = Vre;   Vii = Vri; % mV % initial potential

vpeake = 35e-3; vpeaki = 35e-3; % spike peak amplitude

% numExi=Nei;	numExr=Ner;	numIn=Ni;			numLTS=Nlts;	% num cells to record
total = Ne+Ni;

V=[Vie*ones(Ne,1); Vii*ones(Ni,1);];			% initial v
Vt=[Vte*ones(Ne,1); Vti*ones(Ni,1)];
Vr=[Vre*ones(Ne,1); Vri*ones(Ni,1)];
c=[ce*ones(Ne,1); ci*ones(Ni,1)];
R=[Re*ones(Ne,1); Ri*ones(Ni,1)];
C=[Ce*ones(Ne,1); Ci*ones(Ni,1)];
vpeak=[vpeake*ones(Ne,1) ;  vpeaki*ones(Ni,1)];

% synapses
Sg=zeros(total,total);
Sg(1,2) = 50; Sg(2,1) = 10; % enter as nS

Sg = Sg*1e-9; % convert to S
stdecay = 8; % decay of synapses ms
Stau = ones(size(Sg)).*stdecay; % decay of all synapses ms
SEVr = 0; SIVr = -75e-3; % reversal potential for Ex and In synapses
S_Vr = repmat([SEVr*ones(1,Ne) SIVr*ones(1,Ni)],total,1); % Vr for each synapse (all synapses from a Presynatic cell are in the same column and have the same reversal potential)

% initialize
cg = zeros(size(S)); % current conductance of each synapse

% Setup simulation time
T=200; tau=.25;			% time span and step (ms)
n = round(T/tau);

T1=T/10;

% preallocate intracellular records
vsave = zeros(n,total);
Issave = vsave;

%% run simulation
for i=1:n
    t = i*tau;
    if (t>T1) & (t<T*6/10)
        %         I=-500e-12;
        I=0;
    elseif  (t>T*1/10)
        I=0;
    else
        I=7e-12;
    end;

    if fired
        cg;
    end
    % synpases
    cg(:,fired) = cg(:,fired) + Sg(:,fired); % compute CURRENT synaptic conducance by adding synaptic conductance from fired cells
    Is=sum(cg.*(S_Vr-repmat(V,1,total)),2);				% compute synaptic current from spiked neurons
    cg =  cg.*exp(-tau./Stau); % decay synaptic conductance (exp could be computed outside of loop)

    I=I+Is;


    reset = find(V==vpeak);
    V(reset) = c(reset);% reset

    V = V + tau*(Vr + R.*I - V)./(R.*C*1e3); % euler onestep

    fired=find(V>=Vt & V<vpeak);
    V(fired) = vpeak(fired); % insert spike

    vsave(i,:) = V';
    Issave(i,:) = Is';
end;
% plot(tspan,vsave,[0 T1 T1 max(tspan)],-90+[0 0 10 10]);
figure(99);clf;
subplot(2,1,1);plot([1:n]*tau,vsave); % plot Vm
subplot(2,1,2);plot([1:n]*tau,Issave); % plot synaptic currents

% s = sprintf('
title(num2str(S(1,2)*1e12));
% axis([0 max(tspan) -90 30])