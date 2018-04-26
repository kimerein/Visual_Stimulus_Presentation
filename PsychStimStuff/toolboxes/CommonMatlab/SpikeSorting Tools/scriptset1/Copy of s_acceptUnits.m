% s_acceptUnits.m
% use in s_SpikeSort8...
%     after kmeankluster and plotcluster (or plotmanualcluster)


%Manually enter cluster to accept or join and accept
%%% ACCEPT/JOIN
clear accept_units;
accept_units{1} = {1 4}; %% this number is the cluster or kluster (depending on whic s_plotLusterSpike was last run)
accept_units{2} = {8 19};
accept_units{3} = {14 };
accept_units{4} = {9 23};
accept_units{5} = {10 };
accept_units{6} = {18};  %% removeKluster is not correct if not kluster
accept_units{7} = {12};  %% removeKluster is not correct if not kluster
accept_units{8} = {21};
accept_units{9} = {5};
accept_units{10} = {3 20 7 24};
accept_units{11} = {2};
accept_units{12} = {15};
accept_units{13} = {13};
accept_units{14} = {26};
accept_units{15} = {17};
accept_units{16} = {22};
accept_units{17} = {11 25};
% accept_units{19} = {8};
% accept_units{20} = {27};
% accept_units{20} = {30};
good_type = [1:3 5:17];

% accept_units{10} = {19};
%% REJECT
% reject_unit = [16:20 7 8];
reject_unit =[6 16];

%% Extract Units
unit = struct([]);
for j = 1:size(accept_units,2)
    temp = []; temp1 = temp; temp3 = temp; temp2 = temp;
    for jj = 1:size(accept_units{j},2)
            temp = [temp spikes_to_plot{cell2mat(accept_units{j}(1,jj))}];
            temp1 = [temp1; temp_spikeData{cell2mat(accept_units{j}(1,jj))}];
%             temp2 = [temp2; single_unit{cell2mat(accept_units{j}(1,jj))}];  
%% use if cluster not kluster
            temp2 = [temp2; cell2mat(accept_units{j}(1,jj))];
            if nIntracellsignals >0
                temp3 = [temp3; temp_Intra{cell2mat(accept_units{j}(1,jj))}];
            end
    end
    unit(j).spikes = temp;
    unit(j).spikeData =  temp1;
    unit(j).intraData = temp3';
    if find(good_type == j)
        unit(j).typetag = 1;
    else
        unit(j).typetag = 0;
    end
    %% dont' forget data is sorted by spiketime later
    if j ==1
        removeKluster = temp2;
    else
        removeKluster = [removeKluster; temp2];
    end
end
% removeKluster = [removeKluster; reject_unit']'
% 
% removeSpike = [];keepSpike = [];
% for i = 1: size(klusters,1)
%     if any(removeKluster == klusters(i))
%         removeSpike=[removeSpike i];
%     else
%         keepSpike = [keepSpike i];
%     end
% end
% 
% % recluster w/o
% [PC, SCORE, LATENT, TSQUARE] = princomp(single(allElectrodesXspikes(:,keepSpike)'));
% 
% 

%s_findFirstSpikeInBurst
%% Find spikes that occur less then Xms after another spike of the same
%% unit
X = 50e-3; % 50ms
for i = 1: size(unit,2)
    unit(i).spikeOrder = zeros(1,size(unit(i).spikeData,1));
    temp = unit(i).spikeData(:,1);
  % sort by spiketime
   [temp ix] = sort(temp);
    unit(i).spikes(:,ix(:,1));
    unit(i).spikeData = unit(i).spikeData(ix(:,1),:);
    unit(i).intraData = unit(i).intraData(:,ix(:,1));

    ind_first = find(abs(diff(temp)) > (X*Ts )) + 1 ; 
    unit(i).spikeOrder([1;ind_first]) = 1;
end


dt = output.dt;
save([writedirheader 'Units' '_EE' num2str(triggerEE)],'unit','reject_unit','sourcefile','dt','triggerEE','EEgroup')
