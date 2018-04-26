function createATF(filename, data, colname, varargin)
% function createATF(filename, data, colname, varargin)
%
% creates an atf file from data
% data - each column of data is one sweep
% colname - cell contain a name for each signal
%   since this impementation only support 1 signal and the first
%   col is typically time. These would normally be:
%   colname{1} = 'Time (ms)'; colname{2} = 'Signal (mV)';
% optional input: 'comment','string'
%
% TODO: could be modified to do more than one channel as well as sweeps
% really slow!
% BA 010908
comment = '';

DIR = struct([]);
if nargin>=4
    for i=1:length(varargin)
        if mod(i,2)~=0
            DIR(floor(i/2)+1).param = lower(varargin{i});
        else
            DIR(floor(i/2)).val = varargin{i};
        end
    end
end
for i=1:length(DIR)
    if ~isempty(DIR(i).param)&~isempty(DIR(i).val)
        switch DIR(i).param
            case 'comment'
                comment = DIR(i).val;
            otherwise
        end
    end
end

% check that comment doesn't contain /n
if strfind(comment,'/n'); error('Comment cannot contain newline character "/n"'); end

nsweeps = size(data,2);


%SAVE to atf format
fid = fopen(filename,'w');
fprintf(fid,'ATF 1.0\n1\t%d\n',nsweeps);
fprintf(fid,'"Comment:%s"\n',comment);
fprintf(fid,'"%s"',colname{1}); %%
for j=1:nsweeps-1 %% write colum titles
    fprintf(fid,'\t"%s"',colname{2}); %%
end
fprintf(fid,'\n');
for i=1:size(data,1)
    fprintf(fid,'%f',data(i,1)); %% must not save in %f or %g format so can't use save -ascii
    for j = 2:nsweeps
        fprintf(fid,'\t%d',data(i,j)); %% must not save in %f or %g format so can't use save -ascii
    end
    fprintf(fid,'\n'); %% must not save in %f or %g format so can't use save -ascii
end
fclose(fid)
