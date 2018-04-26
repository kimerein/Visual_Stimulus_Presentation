function imdf = creatdf_f(rawfilename,C,R,F,bsave)
% function imdf = creatdf_f(rawfilename,C,R,F,bsave)
% read in tifs, compute deltaf/f (using average of all frames), scal =1000
% 
% BA 020907

im = loadTILLimg(rawfilename,C,R,F);
%% deltaf 
imdf = deltaf(im,1,F,1000);

if bsave
    %% extract path
    temp = findstr(rawfilename,'\');
    initN = str2num(rawfilename(max(strfind(rawfilename,'_'))+1:max(strfind(rawfilename,'.'))-1));

    writeTILLimg(imdf,[rawfilename(1:max(temp)) 'df_' num2str(initN) '.tif']);
end
