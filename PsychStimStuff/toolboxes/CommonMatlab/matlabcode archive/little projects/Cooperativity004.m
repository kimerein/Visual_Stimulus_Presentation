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
toe = [1:1:5 zeros(1,195)];
Con = Con+ toeplitz(toe);
prob = 0.5; % probability of connection
Thres = [10 20 40 80]; % thresholds
temp = 0;

%close 1
close 3
figure(1)
hold on;
figure(3)
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
forfit = zeros(4,2);
%Watch out i and ii don't match from the loop above to the loop belo
for i = 1 : size(Thres,2) % loop through thresholds
    counter = 0;
    for ii = 1 : NPre % loop through # Total Inputs (
        Prob_Active_Post(i,ii) = sum(Thres(1,i) <  sum(Con(1:ii,:)))/NPost;
        BI_Active_Post(i, ii) = (sum(Prob_Inputs(Thres(1,i):end,ii))); % find all probabilitys where > threshold activations occur
        if (ii > 1 )% && ii < 100*i)
        %    S(i,ii)  = (-BI_Active_Post(i,ii)+ BI_Active_Post(i,ii-1)) *ii/BI_Active_Post(i,ii);
            S(i,ii)  = (-Prob_Active_Post(i,ii)+ Prob_Active_Post(i,ii-1)) *ii/Prob_Active_Post(i,ii);
        end
    end
    %       fitdata(i, 1:2) = polyfit(log10(forfit(:,1,i)),forfit(:,2,i),1)
    x = 1:1:200;
    figure(3)
    plot(x,Prob_Active_Post(i,:),'-b','LineWidth',2);
    plot(x,BI_Active_Post(i,:),'-r','LineWidth',2);

   %  plot(x(2:1:40*i),S(i,1:40*i-1),'-k');
    %plot(x(2:end),S(i,1:size(S,2)-1),'-k');
 


end
