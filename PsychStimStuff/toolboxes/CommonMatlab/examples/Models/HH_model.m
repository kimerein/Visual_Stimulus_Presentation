function yp = HH_model(t,y)
% function dy = HH_model(t,y)
% basic hodgkin huxley model (from sachin code)

global istart istop iappl

if (istart<t)&&(t<istop)
    Iapp = iappl(round(t*10)); %% assume that t is in  ms and iappl contains 10 samples per ms
else
    Iapp = 0; % ?A/cm2
end

% y(1) = V;
% y(2) = h;
% y(3) = n;
% y(4) = a;
% y(5) = b;
yp = zeros(size(y));
% rename to make easy to follow
V = y(1) ;
h = y(2) ;
m = y(3) ;
n = y(4) ;
% 
% C_m = 1;
% Iapp=1; % ?A/cm2
% gNa=215;gK=43;gL=0.813;
% V_Na=215;V_K=-95.0;V_L=-64; %mV

C_m = 1;
gNa=120;gK=36;gL=0.3;
V_Na=45;V_K=-82;V_L=-59; %mV
% # ode's
Vp =	(gNa * m^3 * h * (V_Na - V) + ...
    gK * n^4 * (V_K - V) + gL * (V_L - V) + ...
    Iapp)/ C_m;

 %% not sure why input Vrest (maybe should be reversal?)
% mp=alpham(V,V_L)*(1-m)-betam(V,V_L)*m;
% hp=alphah(V,V_L)*(1-h)-betah(V,V_L)*h;
% np=alphan(V,V_L)*(1-n)-betan(V,V_L)*n;
mp=alpham_HH1(V)*(1-m)-betam_HH1(V)*m;
hp=alphah_HH1(V)*(1-h)-betah_HH1(V)*h;
np=alphan_HH1(V)*(1-n)-betan_HH1(V)*n;

yp(1) = Vp;
yp(2) = hp;
yp(3) = mp;
yp(4) = np;

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