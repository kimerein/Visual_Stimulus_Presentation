x =11
% a = output.data(:,3+x*output.Nchan:3+(x+1)*output.Nchan-2)
% a = output.data(:,3+x*output.Nchan:3+x*output.Nchan+3)
a = output.data(:,7+x*output.Nchan:7+x*output.Nchan+2)
[icasig]= FASTICA(a');
for i = 1:size(icasig,1)
    figure(11);
    subplot(size(icasig,1),1,i)
    plot(-icasig(i,:),'r')
    figure(2);
    subplot(size(icasig,1),1,i)
    plot(a(:,i),'b')
end
icadata = []
for i = 1:output.Nsweeps
    warning  off all
for i = 1:output.Nsweeps
    icadata= [icadata int32(FASTICA(output.data(:,(EEgroup(1)+2)*i*output.Nchan:(EEgroup(1)+2)*i*output.Nchan+length(EEgroup)-1)','verbose','off')'*2^20)];
end
%% indexing wrong getting intracellular data