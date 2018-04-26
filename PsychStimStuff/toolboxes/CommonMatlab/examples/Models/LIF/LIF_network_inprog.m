% inhibitory and excitatory neurons


Ne = 1;         Ni = 1;

Vie = -55e-3;   Vii = -55e-3; % mV % initial potential
Vte=-45e-3;     Vti=-45e-3;	% Threshold
Vre=-44e-3;     Vri=-70e-3; % Resting Potential
ce = -55e-3;    ci = -55e-3; % reset potential
Re = 100e6;     Ri = 20e6; % Ohms 
Ce = 0.1e-9;    Ci = 0.1e-9; %Farads

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

S=zeros(total,total);
S(1,2) = 0; S(2,1) = 0;


T=200; tau=.25;			% time span and step (ms)
n = round(T/tau);

T1=T/10;

% preallocate intracellular records
vsave = zeros(n,total);

%%
for i=1:n
    t = i*tau;
    if (t>T1) & (t<T*6/10)
%         I=-500e-12;
        I=0;
    elseif  (t>T*6/10)
        I=50e-12;
    else
        I=0;
    end;
            Is=sum(S(:,fired),2);				% synaptic current from spiked neurons
            
        I=I+Is;

     %LIF neuron
     reset = find(V==vpeak);
     V(reset) = c(reset);% reset

     V = V + tau*(Vr + R.*I - V)./(R.*C*1e3); % euler onestep
     
     fired=find(V>=Vt & V<vpeak);
     V(fired) = vpeak(fired); % insert spike
     
    vsave(i,:) = V';
end;
% plot(tspan,vsave,[0 T1 T1 max(tspan)],-90+[0 0 10 10]);
plot([1:n]*tau,vsave)
% axis([0 max(tspan) -90 30])