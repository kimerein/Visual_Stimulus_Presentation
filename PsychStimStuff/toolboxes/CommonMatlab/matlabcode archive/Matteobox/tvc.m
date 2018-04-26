function cts = tvc(pars,cps,maxc,logflag) 
% TVC threshold vs contrast function
%
% cts = tvc(pars,cps)
% where cps is a vector of pedestal contrasts, is the threshold expected from the contrast
% response function rvc(pars,cps). 
%
% The threshold is computed numerically, not simply taking the inverse of the derivative.
%
% cts = tvc(pars,cps,maxc) lets you pick if contrasts are bet 0 and 1 (maxc = 1)
% or bet 0 and 100 (maxc = 100, default)
%
% cts = tvc(pars,cps,maxc,'log') takes the log of cts. This is useful for fitting data in log space.
% 
% EXAMPLE:
% cps = logspace(-4,0);
% pars = [20 2 0.2 0.25 1];
% cts = tvc(pars,cps,1);
% loglog(cps,cts);
%
% SKIP TO THE END OF THE FILE TO SEE AN EXAMPLEOF FITTING A TVC
%
% 2000-01 Matteo Carandini
% part of the Matteobox toolbox

if nargin<4
   logflag = 'linear';
end

if nargin<3
   maxc = 100;
end

ccc = maxc*logspace(-4,0,10000);

rrr = rvc(pars, ccc, maxc);

ncps = length(cps);

cts = zeros(1,ncps);
for ic = 1:ncps
   cp = cps(ic);
   icp = find( abs(cp-ccc)==min(abs(cp-ccc)) );
   r2 = rrr(icp)+1;
   if all(rrr<r2)
      ct = maxc; % used to be inf, which is more appropriate but problematic when fitting
   else
      ict = find( abs(rrr-r2)==min(abs(rrr-r2)) );
      ct = ccc(ict)-cp;
   end
   cts(ic) = ct;
end

if strcmp(logflag,'log')
   cts = log10(cts);
end

return

%% ------------------- Code to test the function ----------------------

cps = [0 4 5 8 10 16 20 32 40 64 80];
% Contrasts in percent. Note that there is also a condition cp=0. 
% When plotting in log scale is will be turned into a 0.01

% Measured thresholds
cts = [2.99 0.959 1.63 3.06 6.04 5.04 6.85 7.04 5.27 9.45 9.51];

% guessed parameters
k = 10;
m = 4;
n = 0.5;
sigma = 2;
alpha = 0;

guessedpars = [ k, m, n, sigma, alpha ];

% a bunch of contrasts of the pedestal to draw smooth curves
ccc = logspace(-2,2);

%-- the predicted rvc
rrpred = rvc(guessedpars,ccc);
%-- the predicted tvc
yypred = tvc(guessedpars,ccc)

% when fitting, we'll use the 'log' option, equivalent to calling
% yypred = 10.^tvc(guessedpars,ccc,100,'log')

% Skip to GRAPHICS to give a quick look at how bad the initial guess was

%%	FITS

% let's find the best parameters (forcing alpha to be 0)
% I am not sure alpha should be really used anyway

[err,pars]=fitit('tvc',log10(cts),...
   [1 0.5 0 0.1 0],...
   guessedpars,...
   [100 10 1 20 0],...
	[1 1e-3 1e-3 3 1000 ],cps,100,'log');


%-- the predicted rvc
rrpred = rvc(pars,ccc);
%-- the predicted tvc
yypred = tvc(pars,ccc);

%% 	GRAPHICS

xx = cps;
xx(find(xx==0)) = 0.01;
yypred(find(yypred==100))=NaN; % those were artificially set to 100 by tvc...

figure; ax = [];

ax(1) = subplot(2,1,1)
loglog(xx, cts, 'ko'); hold on
loglog( ccc, yypred, 'k-');
ylabel('Threshold contrast (%)');

ax(2) = subplot(2,1,2);
loglog( ccc, rrpred, 'k-');
ylabel('Response (arbitrary units)');
xlabel( 'Pedestal contrast (%)');

set(ax,'xlim',[0.01 100]);
set(ax,'ylim',[0.01 100]);
set(ax,'dataaspectratio',[1 1 1])


