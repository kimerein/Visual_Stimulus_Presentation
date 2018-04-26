function clusterXT(spikes, useassigns, show, threed, bOFFSETPLOT)
%    clusterXT  temporary utility to show clusters
%      clusterXT(spikes, useassigns, show, threed, bOFFSETPLOT)
%
% Modified by BA

if (nargin < 2 || isempty(useassigns)),  useassigns = spikes.overcluster.assigns;  cmap = spikes.overcluster.colors;  end;
if (nargin < 3 || isempty(show)),  show = unique(useassigns);  end;
if (nargin < 4|| isempty(threed)),  threed = 0; end;

if (nargin < 5|| isempty(bOFFSETPLOT)),  bOFFSETPLOT = 0; end; % BA

if (nargin > 1),
    if (isfield(spikes, 'overcluster') && all(ismember(useassigns, unique(spikes.overcluster.assigns))))
        cmap = spikes.overcluster.colors;
    else
        cmap = jetm(length(show));
    end
end

show(show == 0) = [];

clustlist = unique(useassigns);
t = ((0:size(spikes.waveforms,2)-1)-spikes.threshT)./spikes.Fs;

%%%%%%%%%%%%%%%
cla reset; hold on;
hndl = zeros(1,length(clustlist));
for j = 1:length(clustlist)
    k = clustlist(j);
    members = find(useassigns == k);
    waves = spikes.waveforms(members,:);

	if (k == 0),                color = [0 0 0];
	elseif (ismember(k,show)),  color = cmap(k,:);
    else                        color = Clgy;
	end
	
    if (~isempty(members))
        if (~threed)
           
            if bOFFSETPLOT % BA 
                h = mplot(t, waves+.1*(length(clustlist)-j), 'Color', color);
                text(min(t)*1.1,.1*(length(clustlist)-j),num2str(j)) ;
                set(gca,'YTick',[])
            else
                 h = mplot(t, waves, 'Color', color); % original
            end
            
            set(h, 'ButtonDownFcn', {@raise_me, h});
			if (k == 0), hout = h;
            else         hndl(k) = h;
			end
        else
            [lh,ph] = errorarea(mean(waves,1), std(waves,1,1));
            set(lh, 'Color', brighten(color, -0.6), 'ZData', repmat(k, size(get(lh,'XData'))));
            set(ph, 'FaceColor', color, 'ZData', repmat(k, size(get(ph,'XData'))), 'FaceAlpha', 0.8);
        end
    end
end
hold off; axis tight; xlabel('Time (samples)');  ylabel('Voltage (A/D Levels)');

% if (~threed),  uistack(hndl(show), 'top');
% else  cameratoolbar('SetCoordSys', 'y');
% end

if ((length(show) < 33) && ~threed)
    leg = cell(length(show),1);
    for k = 1:length(show),  leg{k} = num2str(sort(show(k)));  end;
	if (any(useassigns == 0))
		legend([hout, hndl(show)], cat(1, {'Outliers'}, leg), 0);
	else
		legend(hndl(show),leg,0);
	end
end

