function [d dt totalTime] = loadEPfiles(Fdata)
% function [d dt totalTime] = loadEPfiles(Fdata)
%
% loads N fixed episode or gap-free abf files and concatanates them
% INPUT:
%       Fdata struct
%    Fdata(jj).sName  data filename
%             .sweeps Note: (must be EMPTY for gap-free)
%             .start
%             .stop
%             .sChan
%
%  Note: fixed and gap-free cannot be extracted together
%          sampling rate 1/dt must be the same for all files
% BA080507

flagGF = 0; flagFX = 0; d = [];
for jj = 1:length(Fdata)
    fn =Fdata(jj).sName; sweeps = Fdata(jj).sweeps;start = Fdata(jj).start;stop = Fdata(jj).stop;schan = Fdata(jj).sChan;

    if isempty(sweeps)
        sweeps = NaN;
    end
    if isnan(sweeps)%gap free  %% UNTESTED
        if flagFX == 1;
            error('Cannot load fixed abf and gap-free abf at the same time')
        end
        [tmp, si]= EPload(fn,'sweeps',sweeps,'start',start,'stop',stop,'channels',schan);
        d = [d tmp];
        dt = si*1e-6;
        flagGF = 1;
        swT = length(d)*dt ; % s
        Nsw = 1;
    else
        if flagGF == 1;
            error('Cannot load fixed abf and gap-free abf at the same time')
        end
        [tmp, si]= EPload(fn,'sweeps',sweeps,'start',start,'stop',stop,'channels',schan);
        dt = si*1e-6;
        if ischar(stop) | isempty(stop)
            if ~isempty(start);   tmp = tmp(start/dt+1:end,:,:); end;
        else
            tmp = tmp((round(start/dt)+1):round(stop/dt),:,:);
%            tmp = tmp((round(start/dt)+1):round(stop/dt),:,:);
        end
% BA 100608
        if (i>1)&&((size(tmp,2)~=size(d,2) )|| (size(tmp,3)~=size(tmp,3)))
            error('Cannot load abf with different sweep length or channel numbers')
        end
        d = cat(3,d, tmp);

        % get selected time interval (ONLY for fixed episode abf)
        swT = length(d)*dt ; % s
        Nsw = size(d,3);
    end

    if i>1 && dt~=lastdt
        error('Cannot load abfs with different sampling rates')
    end
    lastdt = dt;
end
totalTime= Nsw*swT;