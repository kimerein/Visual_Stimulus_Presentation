% editz.m - A primitive filter design program.
% ---------------------------------------------------------------
% This function demonstrates the usefulness of plt's data editing
% capability. Two plots are created, one showing the poles and
% zeros of a z-plane transfer function and the other showing the
% magnitude and phase of it's frequency response. The frequency
% response plot automatically updates every time you move any of
% the z-plane roots using the mouse or direct keyboard entry.
%
% In the frequency plot, the usual x-cursor edit boxes show the
% cursor location as a fraction of the sample rate. The Xstring
% parameter is used to show this as an angular measure (in degrees)
% just to the right of the x-cursor readout. The usual y-cursor
% edit box shows the magnitude response in dB. The Ystring parameter
% is used to show this in linear form (in percent) just to the
% right of the y-cursor readout. The AxisLink parameter is used so
% that by default the mag/phase axes are controlled separately.
%
% In the pole/zero plot, the usual x and y-cursor edit boxes show
% the pole/zero locations in cartesian form. The Xstring parameter
% is used to show the polar form just to the right of the x-cursor
% readout. Note that the plt command is followed by axis('equal')
% so that the unit circle always looks like a circle even after
% you resize the figure window.

% ----- Author: ----- Paul Mennen
% ----- Email:  ----- paul@mennen.org

function editz()                % A primitive filter design program.
  z = roots([1 4.3 8 8 4.3 1]); % initial numerator polynomial
  p = roots([1 .38 .82]);       % initial denominator polynomial
  z = z(find(imag(z)>-1e-5));   % plot only upper half of unit circle
  p = p(find(imag(p)>-1e-5));
  uc = exp((0:.01:1)*pi*1j);    % uc = half unit circle (101 points)
  x = 0:.001:.5;                % frequency axis
  set(0,'units','pix');  ssz = get(0,'screensize');
  pos = [7 46 755 480];         % figure positions
  if ssz(4) > 1180  pos2 = pos + [0 pos(4)+40 0 0];
  else  pos2 = pos([3 4]);
        pos2 = [ssz([3 4]) - pos2 - pos([1 2]) pos2];
  end;
  % frequency response plot
  S.fr = plt(x,x,x,x,'FigName','Frequency response','Ylim',[-100 1],...
       'TRACEid',{'Mag','Phase'},'LabelY',{'dB' 'Phase'},...
       'LabelX','Fraction of sample rate','Title','QQ','AxisLink',0,...
       'YlimR',[-400 200],'Xstring','sprintf("%4.2f\\circ",360*@XVAL)',...
       'Ystring','sprintf("%6.4f%%",100*min(1,10^(@YVAL/20)))',...
       'ENApre',[0 0],'Position',pos2,'Options','Slider-Y');
  hfig = gcf;
  x1 = -.18;  x2 = -.02;  dx = (x2-x1)/3; % position of spare box
  y1 = -.16;  y2 = -.08;
  zs = complex(x1+dx,(y1+y2)/2);          % append spare zeros and poles
  z = [zs; zs; z];  p = [z(1:2)+dx; p];
  S.pz = plt(z,p,uc,'Xlim',[-1.1 1.1],'Ylim',[-.19,1.3],'ENAcur',[1 1 0],... % pole/zero plot
       'Styles','nn-','Markers','oxn',...
       'LabelX','real','LabelY','imag','FIGname','Pole/Zero plot',...
       'TRACEid',['Zeros '; 'Poles ';'circle'],'Position',pos,...
       'TRACEc',[0 1 1; 1 1 0; .5 .5 .5],'Options','Slider-X-Y',...
       'Xstring','sprintf("  (%4.3f, %4.2f\\circ)",abs(@XY),angle(@XY)*180/pi)');
  axis('equal');
  setappdata(gcf,'DataEdit',0);  % make it easier to edit both x and y coordinates
  set([gcf hfig],'tag','EDITZ','closeReq','delete(findobj(''tag'',''EDITZ''))');
  S.ti = findobj(hfig,'string','QQ');      % find handle of title
  line([x1 x1 x2 x2 x1],[y1 y2 y2 y1 y1]); % draw box for spare poles/zeros
  set(text(x1,y2+.04,' spares'),'color','blue','fontsize',8); % and label it
  S.tx = [...
    text(-.6,1.1, 'Only the roots on or above the x axis are shown','fontweight','bold');
    text(-.6,.70, '- To move a root:');
    text(-.6,.63, '      Click on it');
    text(-.6,.56, '      Then right click on the x-cursor edit box');
    text(-.6,.49, '      Then drag the root to the desired location');
    text(-.6,.42, '      The root will snap to the x-axis if you get close');
    text(-.6,.35, '      A zero will snap to the unit circle if you get close');
    text(-.6,.25, '- To add a root, drag one from the spare box');
    text(-.6,.15, '- To remove a root, drag it below the x-axis')];
  t =  text(-1.46,1,'editz help','user',S.tx,'ButtonDownFcn',@helpTag);
  set([t; S.tx],'units','normal','color',[.7 .7 1]);
  plt('cursor',getappdata(gcf,'cid'),'set','moveCB',{@curCB,S});
  setappdata(gcf,'NewData',0);  curCB(S);  % update the frequency response plot
% end function editz

function helpTag(h,S) % help tag callback
    tx = get(gcbo,'user');
    v = get(tx(1),'visible');
    if v(2)=='n' v='off'; else v='on'; end;
    set(tx,'visible',v); % toggle help visibility
% end function help

function curCB(S) % pz plot cursor callback
  nd = getappdata(gcf,'NewData');
  if isempty(nd) return; end;              % skip if data not modified
  if nd set(S.tx,'visible','off'); end;    % turn off help text after first time
  z = complex(get(S.pz(1),'x'),get(S.pz(1),'y'));  % get zeros and poles from plot
  p = complex(get(S.pz(2),'x'),get(S.pz(2),'y'));
  z3 = z(3:end);  p3 = p(3:end);           % don't use the spares in H(z)
  if z(1) ~= z(2)  z3 = [z(1) z3]; end;    % if a spare was used, create another one
  if p(1) ~= p(2)  p3 = [p(1) p3]; end;
  ytol = diff(get(gca,'ylim'))/80;         % snap to tolerance
  z3(find(imag(z3)<-ytol)) = [];           % delete any roots dragged below the x axis
  p3(find(imag(p3)<-ytol)) = [];
  az = find(abs(imag(z3))<ytol);           % any roots close to the x axis?
  ap = find(abs(imag(p3))<ytol);
  z3(az) = real(z3(az));                   % if so, snap to the x axis
  p3(ap) = real(p3(ap));
  az = find(abs(1-abs(z3))<ytol);          % any zeros close to the unit circle?
  z3(az) = exp(angle(z3(az))*1j);          % if so, snap to the unit circle
  z = [z(2) z(2) z3]; set(S.pz(1),'x',real(z),'y',imag(z)); % in case anything changed
  p = [p(2) p(2) p3]; set(S.pz(2),'x',real(p),'y',imag(p));
  z = [z3 conj(z3(find(imag(z3)>1e-5)))];  % append complex conjugates
  p = [p3 conj(p3(find(imag(p3)>1e-5)))];  % append complex conjugates
  [pn pd] = zp2tf(z.',p.',1);              % compute polynomials
  x = get(S.fr(1),'x');
  H = exp(pi*x*2i);
  pv = polyval(pn,H) ./ polyval(pd,H);     % evaluate polynomials around unit circle
  H = 20 * log10(abs(pv));
  % H=20*log10(abs(freqz(pn,pd,2*pi*x)));  % compute frequency response (db)
  PH = angle(pv)*180/pi;                   % compute frequency response (phase)
  set(S.fr(1),'y',H-max(H));               % update frequency response plot, magnitude
  set(S.fr(2),'y',PH);                     % update frequency response plot, phase
  s = 'H(z) = [';                          % put polynomials in graph title
  p = [pn pd];   l = length(p);  ln = length(pn);
  if l>14 fmt='%4w'; else fmt='%5w'; end;
  for k=1:l  switch k case ln,   t = '] / [';
                      case l,    t = ']';
                      otherwise, t = '  ';
             end;
             s = [s plt('ftoa',fmt,p(k)) t];
  end;
  set(S.ti,'string',s);
  setappdata(gcf,'NewData',[]);            % reset data modified flag
%end function curCB
