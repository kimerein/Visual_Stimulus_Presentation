colororder = ['-.r';'-.g';'-.b';'-.c';'-.m';'-.y'];
global colororder2 =['r-';'g-';'b-';'c-';'m-';'y-';'r-';'g-';'b-';'c-';'m-';'y-';'r-';'g-';'b-';'c-';'m-';'y-';'r-';'g-';'b-';'c-';'m-';'y-'];
colororder3 =['r--';'g--';'b--';'c--';'m--';'y--';'r--';'g--';'b--';'c--';'m--';'y--';'r--';'g--';'b--';'c--';'m--';'y--';'r--';'g--';'b--';'c--';'m--';'y--'];
colorO =['r';'g';'b';'c';'m';'y';'k';'r';'g';'b';'c';'m';'y';'k';];
% colordot =['or';'og';'ob';'oc';'om';'oy';'ok';];
colordot =['.r';'.g';'.b';'.c';'.m';'.y';'.k';'.r';'.g';'.b';'.c';'.m';'.y';'.k';'.r';'.g';'.b';'.c';'.m';'.y';'.k';];
tt = size(colormap,1);
temp = colormap;
my_colors=temp([1:round(tt/10):tt],:);
my_colors = [my_colors; my_colors; my_colors;my_colors;my_colors;my_colors];
my_style =['s','o','x','+','*','.','d','v','^','<','>','p','h'];
% my_syle = [my_style, my_style, my_style, my_style];

% better way of doing it !!!! 
% cmap = jetm(length(NK))