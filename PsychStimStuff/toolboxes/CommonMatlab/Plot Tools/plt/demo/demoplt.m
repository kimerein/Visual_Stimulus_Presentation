% demoplt.m - runs all 22 plt example programs

% ----- Author: ----- Paul Mennen
% ----- Email:  ----- paul@mennen.org

function demoplt()

  fcn = {'plt5'    'plt50'   'pltn'      'pltvbar' 'pltquiv' 'gauss'   ...
         'tasplt'  'trigplt' 'subplt'    'subplt8' 'winplt'  'editz'   ...
         'weight'  'curves'  'pub'       'pltvar'  'pltsq'   'movbar(1)' ...
         'dice(0)' 'bounce'  'circles12' 'wfall'};
  set(0,'units','pixels');  ssz = get(0,'screensize');       % get screen size in pixels
  w = ssz(3);
  figure('menu','none','number','off','share','off',...
         'pos',[w-515 50 510 400],'color',[0,.4,.4],'name','demoplt');
  bx = uicontrol('style','listbox','pos',[5 5 500 292],...
                 'fontsize',10,'foreground',[1 1 0],'background',[0 0 0]);
  x = 10;  y = 400-32;
  for k=1:length(fcn)
     f = fcn{k};
     u = uicontrol('Pos',[x y 55 22],'string',f,'callback',{@OneDemo,f});
     x = x + 63;
     if x > 480  x = 10;  y = y - 32; end;
  end;
  uicontrol('Pos',[x y 118 22],'user',0,'BackG',fliplr(get(u,'BackG')),...
            'string','All Demos','fontsize',8,'callback',{@AllDemos,fcn});
  set(findobj(gcf,'type','uicontrol'),'units','norm');
  if w > 1030  % make figure larger for bigger screens
     set(gcf,'pos',[w-690 50 680 680]);
     set(bx,'fontname','FixedWidth');  % in case Andale Mono doesn't exist 
     set(bx,'fontname','Andale Mono');
  end;
%end function demoplt

function OneDemo(h,arg2,func)
  close(findobj('share','on'));
  set(findobj(gcf,'style','push'),'fontsize',8,'fontw','normal');
  set(findobj(gcf,'string',func),'fontsize',10,'fontw','bold');
  f = fopen(which(func));
  s = {};
  while 1  ln = fgetl(f);
           if ~ischar(ln), break, end
           s = [s; {ln}];
  end
  fclose(f);
  set(findobj(gcf,'style','listbox'),'val',1,'string',s);
  if strcmp(func,'pltvar') evalin('base',func);
  else                     eval(func);
  end;
  set(0,'child',flipud(get(0,'child'))); % force demoplt window on-top
%end function OneDemo

function AllDemos(h,arg2,fcn)
  nf = length(fcn);
  n = get(h,'user') + 1;  % function to start
  if n>nf set(h,'user',0,'string','All Demos');
          close(findobj('share','on'));
          set(findobj(gcf,'style','push'),'fontsize',8,'fontw','normal');
          set(findobj(gcf,'style','listbox'),'string',{'' '' '   --- ALL DEMOS COMPLETED ---'});
  else    set(h,'user',n,'string','click to continue');
          OneDemo(0,0,fcn{n});
  end;
%end function AllDemos
