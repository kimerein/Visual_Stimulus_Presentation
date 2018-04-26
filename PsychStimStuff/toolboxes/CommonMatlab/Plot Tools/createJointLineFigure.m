function createJointLineFigure(A,B,bVertical,fid)
% function createJointLineFigure(A,B,bVertical,fid)
% create plot optimized for figure import into coreldraw
if nargin < 3 ;bVertical = 0; end
if nargin <4
    fid = figure;
else
    figure(fid);
end

if ~bVertical
    set(fid,'Position',[150   150   898/2 937/9]); % right format for figure

    X = [1 2];
    plotJoinLine(A,B,X,3,1);

    plot(B,ones(1,length(B))*X(2),'r.','MarkerSize',35);hold on;
    plot(A,ones(1,length(A))*X(1),'c.','MarkerSize',35);hold on;
    plot(mean(B),X(2),'k.','MarkerSize',35,'linewidth',3)
    plot(mean(A),X(1),'k.','MarkerSize',35,'linewidth',3)
    % xlim([0.5 2.5]);
    set(gca,'YTick',[X])
    % set(gca,'YTickLabel',{'Quartile_2';'Quartile_4'})
    plotset(1)
    % ylim([-.03 1])
    ylim([0.7 2.3])
    % xlabel('psc width (ms)')
else

    set(fid,'Position',[150   150    937/9 898/2]); % right format for figure

    X = [1 2];
    plotJoinLine(A,B,X,3,0);

    plot(ones(1,length(B))*X(2),B,'r.','MarkerSize',35);hold on;
    plot(ones(1,length(A))*X(1),A,'c.','MarkerSize',35);hold on;
    plot(X(2),mean(B),'k.','MarkerSize',35,'linewidth',3)
    plot(X(1),mean(A),'k.','MarkerSize',35,'linewidth',3)
    % xlim([0.5 2.5]);
    set(gca,'XTick',[X])
    % set(gca,'YTickLabel',{'Quartile_2';'Quartile_4'})
    plotset(1)
    % ylim([-.03 1])
    xlim([0.7 2.3])
    % xlabel('psc width (ms)')

end