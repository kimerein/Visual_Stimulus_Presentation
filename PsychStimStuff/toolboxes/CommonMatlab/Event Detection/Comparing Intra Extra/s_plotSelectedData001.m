%% SCRIPT
%% Objective -2
%% loads file intra  extra cellular data from .mat file
%% plots average of selected files from each file on same axis.
%% DEFAULT file selection: file with largest average negative phase of
%%      extracelluar signal
%% Extract data with SpikeExtractionAlignPeaks002.m
colororder2 =['r-';'g-';'b-';'c-';'m-';'k-'];
colororder =['r.';'g.';'b.';'c.';'m.';'k.'];

close all;
% files
readdirheader = 'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Data Analysis\mdata\';
% readdirheader = 'E:\My documents\Scanziani Lab\Data Analysis\mdata\';
readfilename = []; %% for all files

%% all MAT  in directory
path ='*type*.mat';

%GET DATA
[selectDataNum sselectData selectData selectTs] = readMaxSpikeData(readdirheader,path);

proData = baselineCorrect(selectData);
proData = normalize(proData);
%% center data on negative peak
offset = center(proData);
%% align Intra
offset2 = centerIntra(proData);

for(i=1:size(proData,1))
end
%%PLOT
Figure_ID = 300;
timedata = [1:1:size(proData,1)];
index = 0; lastCell = [];
%% Plot selected DATA
tempTs =  max(selectTs(:));  %need because there may be different Ts for different files
nfastb = 0; npyr = 0; spike_fastb=zeros(1,size(proData{1,1},1)*2); spike_pyr = spike_fastb;
for (j = 1:size(sselectData,1) )   %%get selected Data from this cells *.mat file
    %calculate fwhm
    if (~strcmp(lastCell,sselectData{j,1}(strfind(sselectData{j,1}(1,:),'.mat')-3:strfind(sselectData{j,1}(1,:),'.mat')-1)))%% Color each CEll new color
        lastCell = sselectData{j,1}(9:11);
        index = index +1;
    end
    if (index >  size(colororder,1))
        index = 1;
    end

    dfwhm(i) = fwhm(proData{j,1}(:,:,1))/selectTs(j);
    Edfwhm(i) = fwhmPosPhase(proData{j,1}(:,:,3))/selectTs(j);
%     Edfwhm(i) = risetimePosPhase(proData{j,1}(:,:,3))/selectTs(j);
     Ewidth(i) = selectTs(j)*widthnegPhase(proData{j,1}(:,:,3))/selectTs(j);
%     Erisetime(i) = risetimePosPhase(proData{j,1}(:,:,3));  %almost there
%     needs a bit more work
    timedata = [1:size(proData{j,1},1)];

    figure(Figure_ID+100)%intra
    hold on;
    if (strfind(sselectData{j,1},'fastb'))
        scolor2 = '.r';
        scolor = '-r';
        temp = round(size(proData{j,1})/2) -offset(j);
        spike_fastb(temp:size(proData{j,1},1)+temp-1)...
            = spike_fastb(temp:size(proData{j,1},1)+temp-1) +proData{j,1}(:,:,3)';
        nfastb = nfastb +1;
    elseif (strfind(sselectData{j,1},'pyr'))
        scolor2 = '.b'; 
        scolor = '-b';
        spike_pyr(temp:size(proData{j,1},1)+temp-1)...
            = spike_pyr(temp:size(proData{j,1},1)+temp-1) +proData{j,1}(:,:,3)';
        npyr = npyr +1;

    elseif (strfind(sselectData{j,1},'fast'))
        scolor2 = '.k';
        scolor = '-k';
    else
        scolor2 = '.g';
        scolor = '-g';
    end
% scolor2 = '.b';
%     plot(Edfwhm(i) *1000,dfwhm(i)*1000,scolor2,'MarkerSize',20)
    plot(Ewidth(i) *1000,Edfwhm(i)*1000,scolor2,'MarkerSize',20)
    ylabel('FWHM (ms) Extra PostNegPhase PostPhase')
    xlabel('Width (ms)Extra');
    if false
        scolor = colororder2(index,1)
    end
    figure(Figure_ID+200)
        plot(Ewidth(i) *1000,dfwhm(i)*1000,scolor2,'MarkerSize',20)
        hold on;
%     xlabel('FWHM (ms) Extra PostNegPhase PosPhase')
    xlabel('Width (ms) Extra NegPhase')
    ylabel('FWHM (ms)Intra');
    figure(Figure_ID+201)
        plot(Edfwhm(i) *1000,dfwhm(i)*1000,scolor2,'MarkerSize',20)
        hold on;
%     xlabel('FWHM (ms) Extra PostNegPhase PosPhase')
    xlabel('FWHM (ms) Extra PostNegPhase PostPhase')
    ylabel('FWHM (ms)Intra');
if (false)
    figure(Figure_ID)%intra
    %         subplot(2,1,1)
    plot((timedata -offset2(j))/selectTs(j),proData{j,1}(:,1,1),scolor);
    %     line([60,60],[50, -50],'Color','k');
    hold on;

    sSpike= sprintf('%s Sw%d %2.2f E%2.2f',sselectData{j,1},selectDataNum(1,j),dfwhm,Edfwhm);
    stemp2 = sprintf('%s\nMax Negative Phase',path);
    title(stemp2,'Interpreter','none');
%     text(45/tempTs,-j*5,sSpike,'FontSize',8,'Interpreter','none','Color',scolor(2:end));

    figure(Figure_ID+1)%extra
    %         subplot(2,1,1)
    plot((timedata -offset(j))/selectTs(j),proData{j,1}(:,1,3),scolor);
    %     line([60,60],[50, -50],'Color','k');
%     hold on;
%     text(45/tempTs,0-j*.07,sSpike,'FontSize',8,'Interpreter','none','Color',scolor(2:end));
    title(stemp2,'Interpreter','none');
    %     ylim([-1,2])
end
% sSpike
       pause;


end%%END Plot selected DATA

%Align
% temp1 = find(spike_fastb==min(spike_fastb));
% temp2 = find(spike_pyr==min(spike_pyr));
% offset = temp1-temp2;
% temp = circshift(spike_fastb',-offset)/max(spike_fastb)
% temp(end-offset:end) = 0;
% spike_pyr = spike_pyr/max(spike_pyr)
% figure;
% plot(temp)
% hold on;
% plot(spike_pyr,'-r')
% plot(spike_pyr-temp','-g');
% %
