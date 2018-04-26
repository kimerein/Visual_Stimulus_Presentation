function  computeCSD(data,Fs,SITESPACING)
% function  computeCSD(data,Fs,SITESPACING)
%  data - CHN x TIME
%  Fs - sampling rate
% BA 082009
%  SiteSpacing 
% ASSUMES 16 channesl
% see email regarding CSD
% based on 
% Experimental optimization of current source-density technique for anuran cerebellum.
% Freeman JA, Nicholson C.J Neurophysiol. 1975 Mar;38(2):369-82.


bsmooth = 1;

FPOSITION = [2653        -141         520         370];
H = 20003; % figure handle
if nargin < 3
    SITESPACING = 50; %um;
end

try  % if H doesn't exist create it
    get(H,'Name');
catch
    figure(H)
    set(H,'Position',FPOSITION);
%     set(H,'Renderer','OpenGL'); % for some reason this screws up the plot
    set(H,'NumberTitle','off');
    set(H,'Name','CSD');
    lastH = [];
end
set(0,'CurrentFigure',H);



if bsmooth % smooth
    if 1 % smooth with weighted average of nearest neighbors "S1" from Freeman
        % paper
        for i = 2:size(data,1)-1
            smdata(i-1,:) = (data(i-1,:) + data(i+1,:) +2*data(i,:))/4;
        end
%         "D2"
        CSD = -1*(smdata(1:end-2,:) + smdata(3:end,:) - 2*smdata(2:end-1,:))/SITESPACING/SITESPACING;
    else  %% AM NOT SURE BELOW "D4" IS RIGHT!
        % smooth with weighted average of nearest neighbors "S3"
        clear smdata;
        for i = 2:size(data,1)-1;        smdata(i-1,:) = 3*((data(i-1,:) + data(i+1,:)) +4*data(i,:))/10;    end
        for i = 2:size(smdata,1)-1;        sm2data(i-1,:) = 3*((smdata(i-1,:) + smdata(i+1,:)) +4*smdata(i,:))/10;    end % smooth 2x
        smdata = sm2data;
        %"D4"
        CSD = (-5*(smdata(5:end-2,:)+smdata(3:end-4,:))  +6*(smdata(6:end-1,:) + smdata(2:end-5,:)) +9*(smdata(7:end,:) + smdata(1:end-6,:)) - 20*smdata(4:end-3,:))/100/SITESPACING/SITESPACING;
    end
end
% 1 site step size
%         v(x-h) + v(x+h)         - 2 v(x) % see Freeman Nicholson
newc = imresize(CSD,[size(CSD,1)*SITESPACING 5*size(CSD,2)],'bilinear');
T = [1:size(newc,2)]/size(newc,2)*length(data)/Fs; % time
S = [1:size(newc,1)]/size(newc,1)*650 + SITESPACING;  % space
subplot(1,3,1)
imagesc(T,[1:size(newc,1)],newc);
title('1 site step, Smooth nearest neighbor')

% 2 site step size
%         v(x-h) + v(x+h)         - 2 v(x) % see Nicholson
CSD = (data(1:12,:) + data(5:16,:) - 2*data(3:14,:))/(-4)/SITESPACING/SITESPACING; 
newc = imresize(CSD,[size(CSD,1)*SITESPACING 5*size(CSD,2)],'bilinear');
S = [1:size(newc,1)]/size(newc,1)*size(CSD,1)*SITESPACING + SITESPACING;  % space
subplot(1,3,2)
imagesc(T,[1:size(newc,1)],newc);
title('2 site step')

% 3 site step size
%         v(x-h) + v(x+h)         - 2 v(x) % see Nicholson
CSD = (data(1:10,:) + data(7:16,:) - 2*data(4:13,:))/(-9)/SITESPACING/SITESPACING; 
newc = imresize(CSD,[size(CSD,1)*SITESPACING 5*size(CSD,2)],'bilinear');
S = [1:size(newc,1)]/size(newc,1)*size(CSD,1)*SITESPACING + SITESPACING*2;  % space
subplot(1,3,3)
imagesc(T,[1:size(newc,1)],newc);
title('3 site step')

printf('Watch out, DEPTHS may not be right')

% % 1 site step size
% %         v(x-h) + v(x+h)         - 2 v(x) % see Nicholson
% CSD = (data(1:14,:) + data(3:16,:) - 2*data(2:15,:))/SITESPACING/SITESPACING; 
% newc = imresize(CSD,[ 700 5*size(CSD,2)],'bilinear');
% subplot(1,3,1)
% S = [1:size(newc,1)]/size(newc,1)*700;  % space
% imagesc(T,S,newc);
% 
% title('1 site step')
% 

ylabel('depth um')
xlabel('time (s)')
