%% Status: 
% FS works.
% test R and tau
% then run with noise and compare with RS
% (http://senselab.med.yale.edu/ModelDb/showmodel.asp?model=66268&file=\ca1
% _model_db\nca.ode)
% cell or same cell with different resistance and cap


global istart istop iappl
tspan = [0 200]; % time in ms
clear MemV m
istart = 100; istop = 400;
y0 = [-50 .5 .5 .5 .5];
options=odeset('InitialStep',10^(-3),'MaxStep',10^(-1),'stats','on');
nsw = 1;
for sw= 1:nsw
iappl(1:(istop-istop)*10) = -10+3*(sw-1);
[T,MemV] = ode23(@HH_FSmodel,tspan,y0,options);
m{1,sw} = MemV;
m{2,sw} = T;
m{3,sw} = iappl;
end
%% currently FS doesn't seem to do anything
figure(1);clf;
for sw = 1:nsw
    plot(m{2,sw},m{1,sw},'-')
end

%%
options=odeset('InitialStep',10^(-3),'MaxStep',10^(-1),'stats','on');
dt = 1;
tspan = [0 200]; % time in ms
% y0 = [-70 0.5961 0.0529 .3177];
y0 = [-40 0.0303  0.6271 0.7292];
[T,MemV] = ode45(@HH_model,tspan,y0,options);

%% currently FS doesn't seem to do anything
figure(2);clf;plot(T,MemV(:,1),'-')