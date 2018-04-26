function sexptname = createExptname(STF,exptparams)
% function sexptname = createExptname(STF,exptparams)

s_RIGSPECIFIC_SPIKESORT_CONFIG;

temp =dirc(fullfile(DAQSAVEPATH,[STF(1).filename '.*']));
if isempty(temp)
    error('no STF.filename found check path')
end
sexptname = sprintf('%s',datestr(temp{4},29));

for i = 1:length(exptparams) % convert to string
    sexptname = sprintf('%s\t%s',sexptname, num2str(exptparams{i}));
end


