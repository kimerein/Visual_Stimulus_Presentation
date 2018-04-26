colororder = ['-.r';'-.g';'-.b';'-.c';'-.m';'-.y'];
colororder2 =['r-';'g-';'b-';'c-';'m-';'y-';'r-';'g-';'b-';'c-';'m-';'y-';'r-';'g-';'b-';'c-';'m-';'y-';'r-';'g-';'b-';'c-';'m-';'y-'];
colororder3 =['r--';'g--';'b--';'c--';'m--';'y--';'r--';'g--';'b--';'c--';'m--';'y--';'r--';'g--';'b--';'c--';'m--';'y--';'r--';'g--';'b--';'c--';'m--';'y--'];
colorO =['r';'g';'b';'c';'m';'y';'k';];
tt = size(colormap,1);
temp = colormap;
my_colors=temp([1:round(tt/10):tt],:);
my_colors = [my_colors; my_colors; my_colors;my_colors;my_colors;my_colors];
my_style =['.','o','x','+','*','s','d','v','^','<','>','p','h','.','o','x'
    ,'.','o','x','+','*','s','d','v','^','<','>','p','h','.','o','x'];
my_syle = [my_style, my_style, my_style, my_style];