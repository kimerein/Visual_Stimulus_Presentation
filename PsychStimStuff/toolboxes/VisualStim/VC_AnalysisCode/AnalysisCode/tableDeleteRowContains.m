function tableDeleteRowContains(filename,TABLENAME,nCol,valuetoDel)

try
    load(filename,sprintf('%s',TABLENAME));
    % turn everything into a string for comparison (it it was one to begin with it doesn't
    % hurt)
    temp = cellfun(@num2str,eval(sprintf('%s(:,%d)',TABLENAME,nCol)),'UniformOutput',false);
    ind = find(strcmp(temp,valuetoDel));
    if ~isempty(ind)
        eval(sprintf('%s(%d,:)=[]',TABLENAME,ind));
        save(filename,sprintf('%s',TABLENAME),'-append');
        printf('DELETED %d entries TABLE: %s (containing ''%s'')',length(ind),TABLENAME,valuetoDel)     
    end
    
catch ME
    rep = getReport(ME, 'basic')
end