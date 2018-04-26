breport = 1;
bsave = 1;
savetype = 'emf';
exptnum = savefilename;


figure(50)

clf
%% Plot on same graph
bSeperatePlot =0;
if ~bSeperatePlot
    %%% for plotting on seperate graphs
    for i=A:B
        %     for i=1:size(event,2)
        if event(i).V < -70
            IN = 0;
            if i==A
                xcolor = 'r';
            else
                xcolor = 'g';
            end
        else
            IN = 1;
            if i==A
                xcolor = 'c';
            else
                xcolor = 'b';
            end
        end

        set(gcf,'Position',[580 50 1000 700],'PaperUnits','normalized');  orient tall
        subplot(2,2,1)
        [a x] =hist(event(i).peak,100)
        title(['Peak PSC Amplitude' ' ' event(i).sfilename],'Interpreter','none')
        xlabel('pA')
        axis tight
        %         if (~IN)
        %         a = -1*a;
        %     end

        if(~IN)
            stairs(-1*x,a,xcolor,'LineWidth',2);
            xlim([min(-1*x) 1000])
        else
            stairs(x,a,xcolor,'LineWidth',2);
            xlim([min(x) 1000])
        end
        hold all
        subplot(2,2,2)
        [a x] =hist(event(i).isi,100)
        stairs(x,a,xcolor,'LineWidth',2);
        title(['Interevent Interval N=' num2str(length(event(i).isi))] )
        xlabel('ms')
        axis tight
        xlim([0 100])
        hold on

        subplot(2,2,3)
        %% NOT truncating distribution
        [a x] =hist(event(i).RiseTau(find(event(i).RiseTau<100)),200)
        stairs(x,a,xcolor,'LineWidth',2);

        title('PSC \tau_{rise}')
        xlim([0 10])
        xlabel('ms')
        hold on

        subplot(2,2,4)
        [a x] =hist(event(i).DecayTau(find(event(i).DecayTau<100)),200)
        stairs(x,a,xcolor,'LineWidth',2);
        title('PSC \tau_{decay}')
        xlim([0 50])
        xlabel('ms')
        hold on


    end
end
figString = ['Events_' exptnum '_' event(i).sfilename]
if breport
    dirtemp = 'Pairs';
    figdesc = [figString];
    savefigure(writedirheader,dirtemp,figdesc,savetype)
end

bSeperatePlot =0;
if bSeperatePlot
    %%% for plotting on seperate graphs
    for i=1:size(event,2)
        if event(i).V < -70
            IN = 0;
            xcolor = 'r';
        else
            IN = 1;
            xcolor = 'c';
        end

        figure(50+i)
        set(gcf,'Position',[580 50 1000 700],'PaperUnits','normalized');  orient tall
        subplot(2,2,1)
        [a x] =hist(event(i).peak,100)
        title(['Peak PSC Amplitude' ' ' event(i).sfilename ' @' num2str(event(i).V) 'mV'],'Interpreter','none')
        xlabel('pA')
        axis tight
        %         if (~IN)
        %         a = -1*a;
        %     end

        if(~IN)
            bar(-1*x,a,xcolor,'EdgeColor',xcolor);
            xlim([min(-1*x) 1000])
        else
            bar(x,a,xcolor,'EdgeColor',xcolor);
            xlim([min(x) 1000])
        end
        hold all
        subplot(2,2,2)
        [a x] =hist(event(i).isi,100)

        title(['Interevent Interval N=' num2str(length(event(i).isi))] )
        xlabel('ms')
        axis tight
        xlim([0 100])
        hold on

        subplot(2,2,3)
        %% NOT truncating distribution
        [a x] =hist(event(i).RiseTau(find(event(i).RiseTau<100)),200)
        bar(x,a,xcolor,'EdgeColor',xcolor);

        title('PSC \tau_{rise}')
        xlim([0 10])
        xlabel('ms')
        hold on

        subplot(2,2,4)
        [a x] =hist(event(i).DecayTau(find(event(i).DecayTau<100)),200)
        bar(x,a,xcolor,'EdgeColor',xcolor);
        title('PSC \tau_{decay}')
        xlim([0 50])
        xlabel('ms')
        hold on

        figString = ['Events_' exptnum '_' event(i).sfilename '_' num2str(event(i).V) ]
        if breport
            dirtemp = 'Pairs';
            figdesc = [figString];
            savefigure(writedirheader,dirtemp,figdesc,savetype)
        end
    end
end
