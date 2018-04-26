 %%%%%%%%%5
 % Analyze a subset of spikes
%Indices of Spikes to analyze
Filename = [07120502 07130504 07130503 07130501 07150508 07150504 07180503 ...
     07120503 07120501 07130502 07140505 07140504 07140503 07140502 ... 
     07140501 07150510 07150507 07150506 07150505 07150503 07150502 07150501 ...
     07180505 07180504 07180503 07180501 07190504 07190501 07200501]; %Group B
% Filename = [07120503 07120501 07130502 07140505 07140504 07140503 07140502 ... 
%     07140501 07150510 07150507 07150506 07150505 07150503 07150502 07150501 ...
%     07180505 07180504 07180503 07180501 07190504 07190501 07200501]; %Group A
anal_spikesIND = [];
for(i=1:size(Filename,2))
 anal_spikesIND = [anal_spikesIND find(anal_spikes(3,:,2) == Filename(1,i))];
end
% Plot Extracted data
close all;
figure_ind = 50;
colororder =['r.';'g.';'b.';'c.';'m.';'y.';'k.';...
    'ro';'go';'bo';'co';'mo';'yo';'ko';...
    'r.';'g.';'b.';'c.';'m.';'y.';'k.';...
    'ro';'go';'bo';'co';'mo';'yo';'ko';...
    'r.';'g.';'b.';'c.';'m.';'y.';'k.';...
    'ro';'go';'bo';'co';'mo';'yo';'ko';...
    'r.';'g.';'b.';'c.';'m.';'y.';'k.';...
    'ro';'go';'bo';'co';'mo';'yo';'ko';...
    'r.';'g.';'b.';'c.';'m.';'y.';'k.';...
    'ro';'go';'bo';'co';'mo';'yo';'ko'];
Ts = 2.0000e-005;

Z1 = 5;  %Column with Intracellular spike Data 
Z2 = 3;
Y1 = 2;
Y2 = 4;
figure(figure_ind)
subplot(2,1,1)
plot(anal_spikes(:,anal_spikesIND,Z1), colororder(1,:))
hold all;
subplot(2,1,2)
plot(anal_spikes(:,anal_spikesIND,Z2),colororder(1,:))
hold all;


Filename = [415006      420005     5126014     5126017  ]; %LG
anal_spikesIND = [];
for(i=1:size(Filename,2))
 anal_spikesIND = [anal_spikesIND find(anal_spikes(3,:,2) == Filename(1,i))];
end
figure(figure_ind)
subplot(2,1,1)
plot(anal_spikes(:,anal_spikesIND,Z1), colororder(2,:))
hold all;
subplot(2,1,2)
plot(anal_spikes(:,anal_spikesIND,Z2),colororder(2,:))
hold all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
%Subset of timecourse 
% NOTE: this truncates data about spike too (i.e. columns 2 and 4)
% extracted_spikes_OUT = extracted_spikes(1:110,:,:);
% clear extracted_spikes;
% extracted_spikes = extracted_spikes_OUT;
