function ds = fxn_HH1(t,s)

global C gNabar gKbar gLbar ENa EK EL istart istop iappl

ds = zeros(4,1);
ds(1) = alpham_HH1(s(4))*(1-s(1))-betam_HH1(s(4))*s(1);
ds(2) = alphan_HH1(s(4))*(1-s(2))-betan_HH1(s(4))*s(2);
ds(3) = alphah_HH1(s(4))*(1-s(3))-betah_HH1(s(4))*s(3);
gNa = gNabar*(s(1)^3)*s(3);
gK = gKbar*(s(2)^4);
if (istart<t)&&(t<istop)
    ix = iappl;
else
    ix = 0;
end
ds(4) = -(1/C)*(gNa*(s(4)-ENa) + gK*(s(4)-EK) + gLbar*(s(4)-EL)) + ix/C;
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

