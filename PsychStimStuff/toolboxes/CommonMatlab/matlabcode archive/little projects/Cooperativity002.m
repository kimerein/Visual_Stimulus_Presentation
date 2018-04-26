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
Con = Con*5*eye(NPre,NPre))
prob = 0.5; % probability of connection
Thres = [10 20 40 80]; % thresholds
temp = 0;

close 1
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
        Prob_Active_Post(i,ii) = log10(sum(Thres(1,i) <  sum(Con(1:ii,:)))/NPost);
        BI_Active_Post(i, ii) = log10(sum(Prob_Inputs(Thres(1,i):end,ii))); % find all probabilitys where > threshold activations occur
   %Extract Points for fitting
        if ((BI_Active_Post(i, ii) ~= -Inf) && (BI_Active_Post(i, ii) ~= 0) && (counter < 5))  %Use first 4 points for fit
            counter = counter +1;
            forfit(counter,1:2,i) = [ii BI_Active_Post(i, ii)];
       end
    end
 fitdata(i, 1:2) = polyfit(log10(forfit(:,1,i)),forfit(:,2,i),1)
 x = log10(1:1:200);
 figure(3)
   % plot(x,Prob_Active_Post(i,:),'-b','LineWidth',2);    
    plot(x,BI_Active_Post(i,:),'-r','LineWidth',2);
    plot(x,fitdata(i,1)*x+fitdata(i,2),'-k');
  
  
end

