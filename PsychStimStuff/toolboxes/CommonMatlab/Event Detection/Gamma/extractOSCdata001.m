DefineDir
expttype = 'KinateOsc';
writedirpath = [DATAANAL_DIR expttype '\'];
    sumdata = 'SummaryData.mat';
        load([writedirpath sumdata]);
        
        %% find all same cell
        dataindex = [];
        celldate = '2s6c1_2006_07_28';
%         celldate = '2s5c1_2006_07_28';
%         celldate = 's2c1_2006_07_20';
%         date = '2006_07_8';
        for i=1:size(summarydata,2)  %% search for duplicates of the same analysis
            if strfind(summarydata(i).exptnum,celldate)
                dataindex = [dataindex i];
            end
        end

       
     sexpt = [];  clear data;
     figure(99)
     clf
     %% extract selected experiment data
  for i=1:length(dataindex)
      in = dataindex(i);
      Stmean(:,i) = summarydata(in).tmean;
      Semean(:,i) = summarydata(in).emean;
      Sv(i) = summarydata(in).holdVoltage;
      sexpt = [sexpt '\n' summarydata(in).exptnum];
      %% lags
      Slags(i,:) = [summarydata(in).holdVoltage summarydata(in).tmeanlag]
  end
  
  
  %% MAKE PLOTS
  clear data2; clear autoL; clear autoI;
  for i=1:length(dataindex)
      if i ==1
          figure(99);clf; %% unnormalized current
          figure(101);clf; %% xcor
          figure(102);clf; %% powerspectraldensity
          figure(103);clf; %% coherence
      end
      %% SELECT color for holding
      in = dataindex(i);
      if summarydata(in).VC
          if summarydata(in).holdVoltage < -70
              xcolor = '-r';
          elseif summarydata(in).holdVoltage > -20
              xcolor = '-g';
          else
              xcolor = '-c';
          end
      else
          xcolor = '-k';
      end
      
      %% PLOT  %% unnormalized Currents
      figure(99);
      plot([1:size(Stmean,1)]*summarydata(dataindex(i)).dt.*1000,Stmean(:,i)-mean(Stmean(:,i)),xcolor);
      hold on
      axis tight
      xlabel('ms')
      ylabel('pA')

      %% PLOT XCORR
      figure(101)
      data2(:,i) = summarydata(in).xcorrLFPIntra;
      %     autoL(:,i) = summarydata(in).autocorrLFP;
      %     autoI(:,i) = summarydata(in).autocorrIntra;
      plot(([1:size(data2,1)]-size(data2,1)/2 +1).*summarydata(dataindex(i)).dt.*1000,data2(:,i),xcolor);
      hold all

      axis tight; xlabel('ms');
      xlim([-10 10])
      xlim([-100 100])

      %% PLOT %% psd
      figure(102)
      plot(summarydata(in).psdLFP(:,1),summarydata(in).psdLFP(:,2),xcolor)
      xlabel('Hz');  ylabel('arb'); hold on; axis tight;
      xlim([0 100]); 
      
      %% PLOT %% coherence
      figure(103)
      plot(summarydata(in).cohLFPIntra(:,1),summarydata(in).cohLFPIntra(:,2),xcolor)
      xlabel('Hz');  ylabel('arb'); hold on; axis tight;
      xlim([0 100]); 
  end
  
  %%%%%%%%%%%%%%%%%
    %% PLOT conductance
            %% assume we are at the reversal potential, so all the current is Ex/In
            %% ADD CODE to CORRECTION factor if not actually at reversal
    ExRev = 7.5; InRev = -87
  in1(1) = dataindex(find(v<-65,1)); in1(2) = dataindex(find(v> -30,1))
  is(1) = summarydata(in1(1)).holdVoltage- ExRev;
  is(2) = summarydata(in1(2)).holdVoltage - InRev;
  figure(100)
  clf
  for i=1:length(dataindex)
      ind =find(in1== dataindex(i));
      if(ind==1)
          temp = Stmean(:,i)/is(1)
          plot([1:size(Stmean,1)]*summarydata(dataindex(i)).dt.*1000,(temp)-mean(temp),'r')
          hold all
      end
      if(ind==2)
          temp = Stmean(:,i)/is(2)
          plot([1:size(Stmean,1)]*summarydata(dataindex(i)).dt.*1000,temp-mean(temp),'g')
          hold all
      end
  end
  ylabel('\Delta nS') %% because there is ongoing synaptic input and it is hard to say when it is zero.
  xlabel('ms')
  title(['V_Ex :' num2str(ExRev) ' V_In :' num2str(InRev)]);

  axis tight
