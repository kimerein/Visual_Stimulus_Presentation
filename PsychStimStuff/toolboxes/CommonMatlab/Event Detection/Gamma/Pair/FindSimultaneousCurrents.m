w = 6% in ms
ii=0;
out = [];
lastin = 1;
for i = 1: length(ax)
    in = match(ax(i),ay(lastin:min(lastin+20,length(ay))),w) ;

    if in~=-1
        in  = in +lastin -1;
        ii = ii+1;
        out(ii,:) = [ax(i) ay(in)];
        lastin = in;
    end
    
end