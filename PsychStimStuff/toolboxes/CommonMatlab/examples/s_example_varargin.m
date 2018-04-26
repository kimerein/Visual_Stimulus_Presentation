%defaults
srebin = 1;
range =500;

%% example varargin

DIR = struct([]);
if nargin>=6
    for i=1:length(varargin)
        if mod(i,2)~=0
        DIR(floor(i/2)+1).param = varargin{i};
        else 
        DIR(floor(i/2)).val = varargin{i};
        end
    end
end
for i=1:length(DIR)
    if ~isempty(DIR(i).param)&~isempty(DIR(i).val)

        switch DIR(i).param
            case 'range'
                range = DIR(i).val
            case 'bsz'
                srebin = DIR(i).val
            otherwise
        end
    end
end
