%Bassam 7/30/05
%Compute Probability of Postsynpatic neuron firing vs Activated Presynpatic
%neurons.
%  Done 2 ways:
% 1) using a random matrix of connectivity Con.
%    Where rows represent Presynaptic neurons and col Post
% 2) using binomial probablity density function.


%Pre by Post.
NPre = 200;
NPost = 200;
Con = round(rand(NPre,NPost));  % Random Matrix of Connectivity
prob = 0.5; % probability of connection
Thres = [10 20 40 80]; % thresholds
temp = 0;

close 1
close 2
figure(1)
hold on;
figure(2)
hold on;
%Binomial

warning off all
clear Prob_Inputs;
% Find Probability PostSyn cell gets i inputs Given ii TOTAL inputs
% (Binomial)
for ii = 1:NPre % Total Inputs
    for i = 1:ii % Active Inputs
        Prob_Inputs(i,ii) =  nchoosek(ii,i)*prob^(i)*(1-prob)^(ii - i); 
        % Probability Po
    end
end

%Watch out i and ii don't match from the loop above to the loop below
for i = 1 : size(Thres,2) % loop through thresholds

    for ii = 1 : NPre % loop through # Total Inputs (
        Prob_Active_Post(i,ii) = sum(Thres(1,i) <  sum(Con(1:ii,:)))/NPost;
        BI_Active_Post(i, ii) = sum(Prob_Inputs(Thres(1,i):end,ii)); % find all probabilitys where > threshold activations occur
    end
    figure(2)
    loglog(Prob_Active_Post(i,:),'-b');    
    loglog(BI_Active_Post(i,:),'-r');
end


