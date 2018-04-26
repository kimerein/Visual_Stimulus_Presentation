function yp = HH_FSmodel(t,y)
% function dy = HH_FSmodel(t,y)
% Golomb D, Donner K, Shacham L, Shlosberg D, Amitai Y, Hansel D (2007)
% Mechanisms of Firing Patterns in Fast-Spiking Cortical Interneurons. PLoS Comput Biol 3:e156
% y(1) = V;
% y(2) = h;
% y(3) = n;
% y(4) = a;
% y(5) = b;
yp = zeros(size(y));
% rename to make easy to follow
V = y(1) ;
h = y(2) ;
n = y(3) ;
a = y(4) ;
b = y(5);

Iapp=4.35; % ?A/cm2
gA=0.39; % mS/cm2 
theta_m=-24.0; %mV
gNa=112.5;gK=225.0;gL=0.25;
sigma_m=11.5;
theta_h=-58.3  ;sigma_h=-6.7;
theta_n=-12.4  ; sigma_n=6.8;
theta_t_h=-60  ; sigma_t_h=-12.0;
theta_t_na=-14.6; sigma_t_na=-8.6;
theta_t_nb=1.3  ; sigma_t_nb=18.7;
theta_a=-50    ; sigma_a=20;
theta_b=-70    ; sigma_b=-6;
tau_b=150; tau_a=2; % ms?
power_n=2.0;
V_Na=50.0;V_K=-90.0;V_L=-70.0; %mV

m_inf=gammaf(V,theta_m,sigma_m);
h_inf=gammaf(V,theta_h,sigma_h);
n_inf=gammaf(V,theta_n,sigma_n);
a_inf=gammaf(V,theta_a,sigma_a);
b_inf=gammaf(V,theta_b,sigma_b);

% # ode's
Vp=-gNa*m_inf^3*h*(V-V_Na)-gK*(n^power_n)*(V-V_K)-gL*(V-V_L)-gA*a^3*b*(V-V_K)+Iapp;
hp=(h_inf-h)/tau_h(V,theta_t_h,sigma_t_h);
np=(n_inf-n)/tau_n(V,theta_t_na,sigma_t_na,theta_t_nb,sigma_t_nb);
ap=(a_inf-a)/tau_a;
bp=(b_inf-b)/tau_b;

yp(1) = Vp;
yp(2) = hp;
yp(3) = np;
yp(4) = ap;
yp(5) = bp;
