clear all
close all
mystyle =['.','o','x','+','*','s','d','v','^','<','>','p','h','.','o','x'
    ,'.','o','x','+','*','s','d','v','^','<','>','p','h','.','o','x'];
% colororder2 =['r-';'g-';'b-';'c-';'m-';'k-'];
colororder =['r.';'g.';'b.';'c.';'m.';'k.'];
readdirheader = 'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Data Analysis\mdata\';
% readdirheader = 'E:\My documents\Scanziani Lab\Data Analysis\mdata\';

readfilename = []; %% for all files
temp = strcat(readdirheader,'Copy*.mat');
filedata = struct2cell(dir(temp));

if isempty(filedata)
    stemp = sprintf('No Files found at: %s\nCheck path',temp)
    error(stemp);
end

init = 0;
for(j = 1:size(filedata,2))
    j
    readfilename =   filedata{1,j} ;
    fileIN = strcat(readdirheader,readfilename);

    load(fileIN,...
        'processedFiles','Ts','Nspikes_in_File','extracted_spikes', 'misalign_ind', 'EE_coord','-mat');
%     if sum(sum(EE_coord(:,2:end))) ~=0
%         if (size(Nspikes_in_File,2) == size(EE_coord,2))
            % mEE = zeros(size(Nspikes_in_File),size(extracted_spikes,1));
            clear temp11;
            %% there are erros in the .mat files some have more EEcoord
            %% then spikes some more spikes then EEcoord this is not a real
            %% fix..
            temp22 = min(size(EE_coord,1),size(Nspikes_in_File,2)); %% hack see above comment
            for i = 1: temp22
                ii = i + init;
                start_ind = sum(Nspikes_in_File(1,1:i-1));
                % correct for dc offset (bit of a hack)
                mEE(ii,:) = (mean(extracted_spikes(:,start_ind+1:start_ind+Nspikes_in_File(1,i),3),2)'- mean(mean(extracted_spikes(1:20,start_ind+1:start_ind+Nspikes_in_File(1,i),3),1)));
                mIntra(ii,:) = (mean(extracted_spikes(:,start_ind+1:start_ind+Nspikes_in_File(1,i),1),2)');
                distance(ii,1) = sqrt(EE_coord(i,2)^2+ EE_coord(i,3)^2 + EE_coord(i,4)^2);
                temp11(i) = min(mEE(ii,:));
                %     distance(i,1) = sqrt(EE_coord(indMin,2)^2+ EE_coord(indMin,3)^2 + EE_coord(indMin,4)^2);
                %     distance(i,1) = sqrt(EE_coord(indMin,2)^2+ EE_coord(indMin,3)^2 + EE_coord(indMin,4)^2);

            end


            amp(init+1:init+temp22) = temp11./min(temp11);
            figure(1)
            plot(mEE');
            
            figure(2)
            plot(mIntra');

 
            figure(3);
            plot(distance(init+1:init+temp22),amp(init+1:init+temp22),mystyle(j),'MarkerSize',15);
            hold on;
            ylim([0 1.1])

            NpositionPerFile(j) = temp22;
            init = init+size(Nspikes_in_File,2);
            pause;
%         end
%     end
end
% title(fileIN(1+max(strfind(fileIN,'\')):end-4),'Interpreter','none')
% title('O-2_280905_typepyr_B02.mat')
xlabel('Distance um');
ylabel('Amplitude uV');

% normalize to max
figure(99)
y = amp;
x = distance';
p = polyfit(x,y,1);
xf = [0:0.1:max(x)]';
f = polyval(p,xf);;
plot(x,y,'.b',xf,f,'-k')


figure(100)
y = amp;
x = distance';
p = polyfit(x,log(y),1);
xf = [0:0.1:max(x)]';
f = polyval(p,xf);;
plot(x,log(y),'.b',xf,f,'-k')


xx = [1:110];
yy = 0.88*exp(-0.03.*xx); %% exponential fit  % manually calculated from linear fit of log plot
