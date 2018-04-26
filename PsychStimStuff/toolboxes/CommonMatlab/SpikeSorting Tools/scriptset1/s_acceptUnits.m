% s_acceptUnits.m
% use in s_SpikeSort8...
%     after kmeankluster and plotcluster (or plotmanualcluster)

%% Extract Units
unit = struct([]);
for j = 1:size(accept_units,2)
    temp = []; temp1 = temp; temp3 = temp; temp2 = temp; temp4 = temp; 
    for jj = 1:size(accept_units{j},2)
            temp = [temp spikes_to_plot{cell2mat(accept_units{j}(1,jj))}];
            temp1 = [temp1; temp_spikeData{cell2mat(accept_units{j}(1,jj))}];
%             temp2 = [temp2; single_unit{cell2mat(accept_units{j}(1,jj))}];  
%% use if cluster not kluster
            temp2 = [temp2; cell2mat(accept_units{j}(1,jj))];
            if nIntracellsignals >0
                temp3 = [temp3; temp_Intra{cell2mat(accept_units{j}(1,jj))}];
                temp4 = [temp4 tempbINSW{(cell2mat(accept_units{j}(1,jj)))}];
            end
    end
    unit(j).spikes = temp;
    unit(j).spikeData =  temp1;
    unit(j).intraData = temp3';
    unit(j).bIN = temp4;
    if find(good_type == j)
        unit(j).typetag = 1;
    else
        unit(j).typetag = 0;
    end
    %% dont' forget data is sorted by spiketime later
%     if j ==1
%         removeKluster = temp2;
%     else
%         removeKluster = [removeKluster; temp2];
%     end
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
    if nIntracellsignals==0
        unit(i).intraData = -1;
        unit(i).bIN = -1;  % 1 if at inhibitory reversla
    else
        unit(i).intraData = unit(i).intraData(:,ix(:,1));
        unit(i).bIN = unit(i).bIN(1,ix(:,1));  % 1 if at inhibitory reversla
    end

    ind_first = find(abs(diff(temp)) > (X*Ts )) + 1 ;
    unit(i).spikeOrder([1;ind_first]) = 1;
end

if bsave
%     dt = output.dt;
save([writedirheader sunitfilename '_EE' num2str(triggerEE)],'unit','sourcefile','dt','triggerEE','EEgroup')

end