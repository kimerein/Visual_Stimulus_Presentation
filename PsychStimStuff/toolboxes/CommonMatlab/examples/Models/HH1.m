function HH1(vstart0,tmax0,iappl0,istart0,istop0)
% function HH1(vstart0,tmax0,iappl0,istart0,istop0)

% http://www.math.ucdavis.edu/~tlewis/teaching/MAT227/HH1.m
global C gNabar gKbar gLbar ENa EK EL istart istop iappl


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables:
% membrance capacitance  uF/cm^2
C = 1.0; 
% Max Na conductance  mS/cm^2
gNabar = 120;
% Max K conductance  mS/cm^2
gKbar = 36;
% Max leakage conductance  mS/cm^2
gLbar = 0.3;
% Na Reversal Potential  mV
ENa = 45;
% K Reversal Potential  mV
EK = -82;
% Leakage Reversal Potential  mV
EL = -59;

vrest = -70;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% total simulations time  msec
tmax = tmax0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Current pulse parameters  uA/cm^2
iappl = iappl0; 
istart = istart0;
istop = istop0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Values for m, n, and h are set to rest state

v = vrest;

m = alpham_HH1(v)/(alpham_HH1(v) + betam_HH1(v));
n = alphan_HH1(v)/(alphan_HH1(v) + betan_HH1(v));
h = alphah_HH1(v)/(alphah_HH1(v) + betah_HH1(v));

vstart=vstart0;
v=vstart;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the equations using ode45

% Initial Conditions 
s0 = [m n h v];

% Time span
tspan = [0,tmax];

options=odeset('InitialStep',10^(-3),'MaxStep',10^(-1));
[T,S] = ode45(@fxn_HH1,tspan,s0,options);

clf
figure(1)
subplot(2,1,1), plot(T,S(:,4)), xlabel('time'), ylabel('Voltage')
hold on
plot([istart,istop],[-90,-90],'r','linewidth',4)
subplot(2,1,2), plot(T,S(:,[1:3])), legend('m','n','h')



        function a=alpham_HH1(v)
            theta = (v+45)/10;
            if (theta ==0)
                a = 1.0;
            else
                a = 1.0*theta/(1-exp(-theta));
            end
        end
    
        function b = betam_HH1(v)
            b = 4.0*exp(-(v+70)/18);
        end
    
        function a = alphah_HH1(v)
            a = 0.07*exp(-(v+70)/20);
        end
    
        function b = betah_HH1(v)
            b = 1.0/(1+exp(-(v+40)/10));
        end
    
        function a=alphan_HH1(v)
            theta = (v+60)/10;
            if (theta ==0)
                a = 0.1;
            else
                a = 0.1*theta/(1-exp(-theta));
            end
        end

        function b = betan_HH1(v)
            b = 0.125*exp(-(v+70)/80);
        end
    
        

    
    
end