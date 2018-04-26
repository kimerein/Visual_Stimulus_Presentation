%% Calculate how many cells achieve threshold when they recieve input from a set of Fibers
%% A and B combined, or just A.  
%%  Results in "out"
%%  
%% Amplitude pdf is used to creat Bag that amplitudes are drawn from 

%% creat distribution of amplitudes;
AmpThres = 280; %% threshold
p = [1,20,1.85]; % Amplitude distribution parameters
x = [1:1:300];
tball = 10000; % total "balls in bag"
P = p(1).*exp(-(log(x./p(2))./p(3)).^2); % probabilty of each size ball
% P = p(1).*exp(-(log(x./p(2)).^2./p(3))); % probabilty of each size ball
P = P ./sum(P); %normalize
Bag = [];
for i=1:size(P,2)
    Bag = [Bag i.*ones(1,round(P(i)*(tball-1)+1))];
end

%% Anatomy Paramters
NF = 10000; % number of Fibre
NC = 1000;  % number of Cells
m = 10; % mean number of synapes from Fibres on to cells
s = 1;% variance in connetivity.  More variance more divergence

Fact  = 1;% number of fibers activated by A or B


C = zeros(NC,NF,'single'); % connectivity matrix
out = zeros(NC/5,5);
ind = 0;
for m=1:10:NC/3  % for different cells/fiber
%     tic
% m =1000;  % for sanity check
    NS = m; % draw number of synapse from normal dist
    try
        % find NS connected cells for each of the NF fibers
        cells = reshape(int32((rand(m,NF)*(NC-1)+1) + ones(m,1,'single')*single([0:NF-1].*NC)),1,m*NF);
        C(cells) = Bag(round(rand(size(cells)).*(size(Bag,2)-1)+1)); % pick amplitude
    catch
        i
    end
    %     toc

    clear frac;clear active; clear activeboth;
    %     tic
    for i = 1:10
        % pick fibers activated by A and B
        A = []; B=A;
        while length(unique(A)) < Fact
            A = round(rand(Fact,1)*(NF-1)) +1;
        end
        while length(unique(B)) < Fact
            B = round(rand(Fact,1)*(NF-1))+1;
        end
        frac(i) = sum(C(:,A)&C(:,B))/(sum(C(:,A)) + sum(C(:,B)));
        activeboth(i) = sum(sum(C(:,A)+C(:,B),2)>AmpThres);% above threshold
        active(i) = sum(sum(C(:,A),2)>AmpThres);% above threshold
    end
    % toc
    if mod(m,33) == 0
        [m s mean(frac) mean(active) mean(activeboth)]
    end
    ind = ind +1;
    out(ind,:) = [m s mean(frac) mean(active) mean(activeboth)];
end

out = out(1:find(out(:,1)==0,1,'first')-1,:);
figure(11)
hold off
plot(out(:,5),'r')
hold on
plot(out(:,4),'b')
title(['Threshold =' num2str(AmpThres)]);
xlabel('Cells/Fiber')
ylabel('frac Cells > Threshold')