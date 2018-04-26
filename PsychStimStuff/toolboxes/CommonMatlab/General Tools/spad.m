function str =  spad(str,N)
% function str =  spad(str,N)
 %% pad string with spaces to be N in length
 l = N-length(str);
 if l<0
 warning('str is longer then N')
 elseif l>0
    str = [str  repmat(' ',1,l)];
 end