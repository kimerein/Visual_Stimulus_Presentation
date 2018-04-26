   stemp = 'dataset('
    for i = 1:size(tempcell,2)
        if i>1
            stemp = [stemp ','];
        end
        stemp = sprintf('%stempcell(2,%d)',stemp,i);
    end
    stemp = [stemp ')'];
    tempdataset =  eval(stemp)
    
        
    %     %     make string to eval to create dataset (can't figure out how to pass
    %     %     in an array of both varnames and values (Really fucking frustrating!)
    %     stemp = 'dataset({tempcell(2,:)'
    %     for i = 1:size(tempcell,2)
    %         stemp = sprintf('%s,tempcell{1,%d}',stemp,i)
    %     end
    %     stemp = [stemp '})'];
    %     tempdataset =  eval(stemp)
    %
    
 