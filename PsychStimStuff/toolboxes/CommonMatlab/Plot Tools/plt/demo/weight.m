% weight.m ----------------------------------------------------------
% - The SubPlot argument is used to create an upper and lower axis.
%   The lower axis contains three traces showing the magnitude of 3
%   different weighting functions used in sound level meters (as
%   defined by IEC 651). The upper axis also contains 3 traces showing
%   the same 3 functions except in dB instead of the linear units used
%   for the lower axis.
% - Normally plt only puts one trace on each subplot except for the
%   main (lower) axis. So in this case (with 6 traces) plt puts 5
%   traces on the lower axis and one on the upper. Since we really
%   want 3 and 3, a separate command is used to move two of the traces
%   from the lower axis to the upper.
% - Note that when you cursor any trace in the lower axis, you get
%   the y-axis readout in both linear and dB units, and the cursor
%   in the upper axis automatically moves to the same trace and the
%   same x position. The subplot argument usually will cause plt to
%   do this synchronization, however since that only deals with a single
%   trace for the sublots we have to use the 'moveCB' cursor callback
%   (and case 1 of the switch command) to do this synchronization.
% - The traceID callback ('TIDcback') insures that the traceID box
%   controls the upper axis traces as well (case 2 of the switch)
% - Note the LineWidth argument in the plt call. This illustrates how
%   any line property may be included in the calling sequence. 


% ----- Author: ----- Paul Mennen
% ----- Email:  ----- paul@mennen.org

function weight()
  f = 1000*10.^((-20:13)./10);                    % x axis for plot
  g = f.^2;
  Wc = g ./ ((g + 20.6^2) .* (g + 12200^2));             % C weight
  Wb = Wc .* sqrt(g ./ (g + 158.5^2));                   % B weight
  Wa = Wc .* g ./ sqrt((g + 107.7^2) .* (g + 737.9^2));  % A weight
  ref = find(f==1000);      % use 1000 Hz as the reference frequency
  W = [Wa/Wa(ref); Wb/Wb(ref); Wc/Wc(ref)]; % normalize to reference
  ctrace = [0 1 0; 1 0 1; 0 1 1]; % color for A/B/C weighting respectively
  S.lh = plt(f,[W; 20*log10(W)],'FigName','A B & C Weighting','EnaPre',[0 0],...
       'LineWidth',{2 2 2 1 1 1},'TIDcback',@tidCB,'SubPlot',[50 50],...
       'ENAcur',[1 1 1 0 0],'Xlim',f([1 end]),'Ylim',[-.05 1.25],...
       'TRACEid',{'A weight' 'B weight' 'C weight' ' ' ' '},...
       'Options','Xlog-X-Y','LabelX','Hz','LabelY',{'Magnitude' 'dB'},...
       'FigBKc',[0 0 .15],'PltBKc',[.15 0 0],...
       'Ctrace',[ctrace;ctrace],'Position',[30 50 900 650]);
  S.ax  = getappdata(gcf,'axis');
  S.cid = getappdata(gcf,'cid');
  set(S.ax(2),'ylim',[-60 5]);  plt('grid',S.ax(2));
  set(S.lh(4:5),'parent',S.ax(2));
  plt('cursor',S.cid(1),'set','moveCB',{@curCB,S});
  set(gcf,'user',S.lh);  % save line handles for traceID callback

function curCB(S) % cursor callback (moveCB)
  xy  = plt('cursor',S.cid(1),'get','position');
  ydB = 20*log10(imag(xy));
  [n h] = plt('cursor',S.cid(1),'get','activeLine');
  c = get(h,'color');
  e = findobj(gcf,'style','edit');
  set(e(2),'str',sprintf('%4.2fdB',ydB),'Backgr',c)
  set(findobj(S.ax(2),'markersize',8),'x',real(xy),'y',ydB,'color',c);
%end function curCB


function tidCB()    % traceID callback
  lh = get(gcf,'user');
  set(lh(4:6),{'vis'},get(lh(1:3),{'vis'}));
% end function tidCB
