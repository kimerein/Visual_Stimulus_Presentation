
tau = 0.25; tspan = 0:tau:100;
T1=tspan(end)/10;
V = -55e-3;
Vt=-45e-3;	Vr=-44e-3; c = -55e-3; % reset potential
vpeak = 35e-3; 
vsave = zeros(1,length(tspan));
R = 100e6; C = 0.3e-9;
vsave = [];
for t=tspan
    if (t>T1) & (t<tspan(end)*6/10)
        I=-50e-12;
    elseif  (t>tspan(end)*6/10)
        I=100e-12;
    else
        I=0;
    end;
%     V = V + tau*(0.04*V^2+5*V+140-u+I)/1;
     %LIF neuron
     V = V + tau*(Vr + R.*I - V)./(R.*C*1e3); % euler onestep
     
     fired=find(V>=Vt & V<vpeak);
     
     reset = find(V==vpeak);
     V(reset) = c;% reset
     V(fired) = vpeak; % insert spike
     
    vsave(end+1) = V;
end;
% plot(tspan,vsave,[0 T1 T1 max(tspan)],-90+[0 0 10 10]);
plot(tspan,vsave)
% axis([0 max(tspan) -90 30])