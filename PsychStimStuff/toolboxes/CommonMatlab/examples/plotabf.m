function []= plotabf(readfileheader,readfilenumber,sweep,Chan,varargin)
% function plots Axon abf files
% readfileheader is a string
% readfilenumber is a vector of filenumbers
%% default Parameters
% Will extract all sweeps: sweep = -1 
if nargin <2 || nargin >4
        error('2 to 4 input arguments are required');
end
if nargin <3
    sweep = -1;
end
if nargin < 4
        Chan = [1 3]; %Channels to extract
end
Ts = 1/0.02e-3 ;%Sample rate 50kHz
readdirheader =  'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\'; 
figure_ID = 10;

close all;
for i = 1: size(readfilenumber,1)
    [filename] = createaxonfilename(readfileheader,readfilenumber(i));
    pathfilename = strcat(readdirheader,filename);
    clear my_Data;
    clear N_sweeps;
    [my_Data N_sweeps] = import_abf(pathfilename,sweep,1/Ts,Chan);
    figure(figure_ID+i);
    for j = 1:size(Chan,2)  %Plot Select Channels
%         subplot(2,ceil(size(Chan,2)/2),j);
        plot(my_Data(:,j+1:size(Chan,2):end),'r-');
        hold on;
        title(filename,'Interpreter','none');
%         if (j-1)/2 ==0 %IntraAxis
%             ylim([-100, 60])
%         else %ExtraAxis
%             ylim([-20,20])
%         end
% 
    end
end