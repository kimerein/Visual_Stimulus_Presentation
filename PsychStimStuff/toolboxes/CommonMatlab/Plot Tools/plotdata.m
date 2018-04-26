function hfid = plotdata(data,dt,varargin)
% function hfid = plotdata(data,dt,varargin)
% plots data against time
% INPUT: 
%        data = array of data: each column should be its own channel
%        dt =  time step 
%        'trange', [mintime maxtime] to plot <all data> (in sec)
%        'chn' , chn = vector of channels to plot
%        'sep', scalar (0/1)  <1> plots each channel it its own (seperate) figure.
%        'fid', scalar handle of first figure <1>
%        'color', scalar handle of first figure <1>
%
% BA012706
fid = 1;
bsep = 0;
if isrowvector(data)
    data = data';
end
chn = [1:size(data,2)];
trange = [1*dt size(data,1)*dt];

scolor = '-';
DIR = struct([]);
if nargin>=3
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
            case 'trange'
                trange = DIR(i).val;
            case 'chn'
                chn = DIR(i).val;
            case 'sep'
                bsep = DIR(i).val;
            case 'fid'
                fid = DIR(i).val;
            case 'color'
                scolor = DIR(i).val;
            otherwise
        end
    end
end
% make row vector if col vector

hfid = figure(fid);
for i=1:size(data,2)
    if any(find(chn ==i))

        if bsep
            hfid = figure(fid);
            fid = fid+1;
        end
        ind = int32(floor(trange./dt));
        if ind(2)> size(data,1)
            ind(2) = size(data,1);
        end
        if ind(1)<= 0 
            ind(1) =1;
        end
        plot([double(ind(1)):double(ind(2))].*dt.*1000, data(ind(1):ind(2),i),scolor)
        hold all
        xlabel('time (ms)')
        %% ADD legend
    end
end
