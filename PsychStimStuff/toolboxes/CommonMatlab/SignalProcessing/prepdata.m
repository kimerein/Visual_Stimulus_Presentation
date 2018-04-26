function df = prepdata (data,bplot,dt,dfilter)
% function df = prepdata (data,bplot,dt,dfilter) OR
% function df = prepdata (data,bplot,dt,dfilter (structure) )
%
% Preprocess data (detrend) >200
% so that there is little drift in baseline to allow effective thresholding
% of currents
%
% INPUT :
%          data: each ROW should contain one sweep or channel
%          dfilter =[HIGH_PASS_Hz LOW_PASS_Hz  60Hz_removal bdetrend]
%                   e.g. [20 200 0 0] bandpass 20-200 and no 60Hz removal,
%                   no detrend
%                        [NaN NaN 1]  finds significant peaks (determined by hard-coded p-value) and subtracts
%                        them, detrend
%                        [NaN NaN 40 1] removes significant peak at 40Hz
%                        (note: for all  60Hz_removal values >0 and not
%                        equal to 1,  signifiant power at THAT frequency
%                        is removed, detrend
%   OR 
%         dfilter.filttype = 1 for filtering with butterworth
%                          = 2 for fft (Throw away unwanted frequencies)
% OUTPUT:
%         df = N by M row. Where N is the sweep/channels and M is the
%         timeseries. i.e. each row is a timeseries
% BA061407
sp = 'y';

if nargin < 3
    dt = 1;
end
if nargin < 2
    bplot = 1;
end
if bplot ==0
    sp = 'n';
end
if ~isrowvector(data) % catches the case of column vector having been input instead or row vector
    if size(data,2) ==1
        data = data';
    end
end
if ~isstruct(dfilter)
if length(dfilter)< 4 % for backwards compatibility
    dfilter(4) = 1;
end
filttype = 1; % use filterdata function
else
    if isfield(dfilter,'filttype')
        filttype = dfilter.filttype;
    else
        filttype = 1;
    end
    dfilter = dfilter.dfilter;
end
    
h = NaN;k=1;
for i = 1: size(data,1)
    if dfilter(4)
        df(i,:) = detrend(data(i,:));
    else
        df(i,:) = data(i,:);
    end
    if dfilter(3)
        params.Fs = 1/dt;
        params.fpass =  [1 100];
        params.tapers = [3 5];
        params.pad = 2;
        %         params.err = [1 fs/0.500];
        %         params.f0 = [];
        p = 5/length(data);
        %     p = [];
        %         f0 = [60]; %[60 88]
        if dfilter(3) == 1; f0 = []; else f0 = dfilter(3); end
        
        if ~isempty(f0);
            for j = 1:length(f0) % allows for multiple f0s but there is NO way to specify them when calling the function
                [df(i,:) h(k)] =     rmlinesc(df(i,:),params,p,sp,f0(j));
                k =k+1;
            end
        else
            k =1;
            [df(i,:) h(k)] =     rmlinesc(df(i,:),params,p,sp,f0);
        end
    end
    % TODO should make into bandpass instead of low and hig pass
    switch filttype
        case 1
            if ~isnan(dfilter(1))&(dfilter(1)>0)
                df(i,:) = filterdata(df(i,:),dt,dfilter(1),1);
            end
            if ~isnan(dfilter(2))
                df(i,:) = filterdata(df(i,:),dt,dfilter(2),0);
            end
        case 2 % use fftfilter
            if ~isnan(dfilter(1))&(dfilter(1)>0)
                df(i,:)= fftFilter(df(i,:)',1/dt,dfilter(1),2)';
            end
            if ~isnan(dfilter(2))
                df(i,:)= fftFilter(df(i,:)',1/dt,dfilter(2),1)';
            end          
    end
end
if bplot
    t = 5;% time in sec
    sw = 1; % if more than 1 sw/episode only plot first
    h=figure;
    clf
    plotdata(detrend(data(sw,:),'constant')',dt,'trange',[0 t],'fid',h);
    hold all
    plotdata(df(sw,:)',dt,'trange',[0 t],'fid',h);
    set(h,'Position',[7   538   560   420]);
    % plotdata(df+std(df)*5,dt,'trange',[0 t],'fid',1)
else
    if ~isnan(h); close(h); end
end
