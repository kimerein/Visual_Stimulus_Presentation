% or MuPAD engine
evalin(symengine,'simplify(sqrt(n*p*(1-p))/n*p)')
evalin(symengine,'evalf(int(1/sqrt(2*pi)*exp(-x^2/2),x=-inf..2))')

maple('solve(F=1/(s*A)*exp((V-Vm)^2/2/s^2),s)')
A = sqrt(2*pi); Vm = -55; V = -46; F = 1/NPYR;
maple('NPYR := 1e5;A := sqrt(2*pi); Vm := -55; V := -45; F := 1/NPYR; SS:=evalf(solve(F=1/(S*A)*exp(-(V-Vm)^2/2/S^2),S))')

sE = 'V = Vr-(Vr-Vmax)*exp(-t)';
sD = sprintf('V:=%1.2f;Vr:= %1.2f;Vmax:= %1.2f%;',VTHRES,Vrest,VPYRmax);
maple(sD);
t=str2num(maple(['solve(' sE ',t)']))

% % Checking calculations are correct
% maple('restart;')
% sD = sprintf('t:=%1.2f;Vr:= %1.2f;Vmax:= %1.2f%;T:=1',t,Vrest,VPYRmax);
% maple(sD);
% maple(['solve(' sE ',V)'])
% 
% maple('restart;')
% sD = sprintf('t:=%1.2f;Vr:= %1.2f;V:= %1.2f%;T:=1',t,Vrest,VTHRES);
% maple(sD);
% maple(['solve(' sE ',Vmax)'])


maple('restart;')
sD = sprintf('t:=%1.2f;Vr:= %1.2f;V:= %1.2f%;T:=1',t+tWin,Vrest,VTHRES);
maple(sD);
Vmin = str2num(maple(['solve(' sE ',Vmax)']))


----------------1/(s*A)*exp((V-Vm)^2/2/s^2),s)

1/2/sqrt(2*pi)*exp(-(-46+55))^2/2/4


sD = 'S:=((1-p)^(n-i))/((1-p)^(2*n-i))';
maple(sD)
SD =  'A :=(n!/(k!*(n-k)!))/((2*n)!/(k!*(2*n-k)!))';
maple(SD)
maple('simplify(A)')

str2num(maple(['simplify(' sD ')']))