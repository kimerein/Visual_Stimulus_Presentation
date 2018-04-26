
tau = 0.25; tspan = 0:tau:100;
T1=tspan(end)/10;
V = -55e-3; % mV
Vt=-45e-3;	Vr=-44e-3; c = -55e-3; % reset potential
vpeak = 35e-3; 
vsave = zeros(1,length(tspan));
R = 100e6; C = 0.1e-9; % Ohms and Farads
vsave = [];
for t=tspan
    if (t>T1) & (t<tspan(end)*6/10)
        I=-500e-12;
    elseif  (t>tspan(end)*6/10)
        I=200e-12;
    else
        I=0;
    end;
     %LIF neuron
     reset = find(V==vpeak);
     V(reset) = c;% reset

     V = V + tau*(Vr + R.*I - V)./(R.*C*1e3); % euler onestep
     
     fired=find(V>=Vt & V<vpeak);
     V(fired) = vpeak; % insert spike
     
    vsave(end+1) = V;
end;
% plot(tspan,vsave,[0 T1 T1 max(tspan)],-90+[0 0 10 10]);
plot(tspan,vsave)
% axis([0 max(tspan) -90 30])