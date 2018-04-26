load Unit file
%% fix old unit file
EEgroup = [6 7 8]
for i = 1:size(unit,2)
    unit(i).bIN = ones(1,size(unit(i).spikes,2))
end
triggerEE = EEgroup(end);
save([writedirheader 'Units' '_EE' num2str(triggerEE)],'unit','reject_unit','sourcefile','dt','triggerEE','EEgroup');
%%%%%

defineDir
bfirst =1;
Ts = 1/dt;
breport = 1;
writedirheader = [DATAANAL_DIR sourcefile(1:end-4) '\' 'EE' strrep(num2str(EEgroup),'  ','_') '\'];

% selectUnits = [1:3];
selectUnits = [1:size(unit,2)];
sunitfilename = 'Unit';
s_MSPairFigure
