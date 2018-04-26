%% Calculate fraction (of all cells activated by set of Fibers A & B) which overlap

NF = 10000; % number of Fibre
NC = 1000;  % number of Cels
m = 10; % mean number of synapes from Fibres on to cells
s = 5;% variance in connetivity.  More variance more divergence
C = zeros(NC,NF,'int16'); % connectivity matrix
out = zeros(NC/5,2);
ind = 0;
for s=1:5:NC/5  % for different variances
    for i=1:NF
        NS = round(normrnd(m,s)); % draw number of synapse from normal dist
        if NS> NC;
            NS = NC
        end
        try
        C(round(rand(NS,1)*(NC-1)+1),i) = 1; % pick cells that are contacted to fiber i
        catch
            i
        end
    end

    Fact  = 10;% number of fibers activated by A or B
    clear frac;
    for i = 1:10
        % pick fibers activated by A and B
        A = round(rand(Fact,1)*(NF-1)) +1;
        B = round(rand(Fact,1)*(NF-1))+1;
        frac(i) = sum(C(:,A)&C(:,B))/(sum(C(:,A)) + sum(C(:,B)));
    end

    if mod(s,10) == 0
        [s mean(frac)]
    end
    ind = ind +1;
    out(ind,:) = [s mean(frac)];
end