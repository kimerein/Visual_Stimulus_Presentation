%intracellular movie
if true
    subplot(1,1,1);
    scaleI=10;
    scaleU=25;
    cells=neighbors(2); %total-2 ;%extneigh(1);
    time=1:ceil(((i-1)/tausPerRecSamp));

    %cells are rows, time is columns
    samps=ceil(((i-1)/tausPerRecSamp));
    recording=[exi(:,1:samps);exr(:,1:samps);in(:,1:samps);lts(:,1:samps)];
    scale=300;
    frame=1;
    for k=1:15:size(time,2)
        plot(time(1:k).*tau,lts(extneigh(20),1:k),'r');hold on;title(stemp);
        F(frame)=getframe;
        frame=frame+1;
        xlim([0 max(time) * tau]);
        ylim([-50 50]);
    end;
%     movie(F);
%     movie2avi(F,'spiketrain.avi');
end;


 
