%% documented in:
%%                   Discriminate FS from PYR.cdr
%%
%% originally used for McKnight grant (4/20/05)
%%
%% data read in with "plotSelectedData2.m"
% then run this first for loop
for i = 1:size(proData,1)
    a = proData{i,1}(:,1,3);
    FSparam(i,:) = xSpikeParams(a');
end

%% Extract spikes
%% Extract Parameters
%% data read in with "plotSelectedData2.m"
for i = 1:size(proData,1)
    a = proData{i,1}(:,1,3);
    PYparam(i,:) = xSpikeParams(a');
 end
%        
%         [ 1  2  3  4   5   6   7   8   9  10]
%% parm = [P1 T1 P2 P1W T1W P2W T1F T1R P2R P2F];
figure(106)
%           plot((FSparam(:,8)+FSparam(:,6)).*.002.*abs(FSparam(:,2)./FSparam(:,3)),1,'.b',(PYparam(:,8)+PYparam(:,6)).*.002.*abs(PYparam(:,2)./PYparam(:,3)),1,'.r')  %102
%           plot((FSparam(:,8)+FSparam(:,6)).*.002,0,'.b',(PYparam(:,8)+PYparam(:,6)).*.002,0,'.r','MarkerSize',20)  %102
plot((FSparam(:,5)).*.002,0,'.b',(PYparam(:,10)).*.002,0,'.r','MarkerSize',20)  %102
ylim([-0.2 0.2])

params = [FSparam;PYparam];
figure(109)
plot((params(:,8)).*.002,zeros(size(params,1),1),'.b','MarkerSize',20)  %102
