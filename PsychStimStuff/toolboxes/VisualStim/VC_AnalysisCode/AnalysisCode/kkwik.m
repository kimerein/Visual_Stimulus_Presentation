function [assign nk]= kkwik(data,Fname)
% function kkwik(data,filename)
% run klustakwik on data
%       generates feature (.fet) file with data provided and runs
%       klustakwik on it.
% INPUT: 
%         data - N x M matrix  where N is number of instances and M number
%         of dimensions
%         Fname - path and filename where data should be saved.
% OUTPUT:
%        assign - number of kluster for each instants (vector of length N)
%        ndim - number of dimensions
% BA052009

% ***  DOCUMENTATION from KlustaKwik readme
% The feature file should have a name like FILE.fet.n, where FILE is any string,
% and n is a number.  The program is invoked by running "KlustaKwik FILE n", and
% will create a cluster file FILE.clu.n and a log file FILE.klg.n.  The number n
% doesn't serve any purpose other than to let you have several files with the
% same file base.
% 
% The first line of the feature file should be the number of input dimensions. 
% The following lines are the data, with each line being one data instance,
% consisting of a list of numbers separated by spaces.  An example file test.fet.1 is provided. Please note that the features can be the sample values of a putative waveform event.


KKPATH = 'C:\KlustaKwik\KlustaKwik-1_7.exe';
temp = strfind(Fname,'\');
SAVEPATH = Fname(1:temp(end));
% ADD deal with multiple .n
SAVEFILENAME = [Fname(temp(end)+1:end) '.fet.1'];

fid = fopen([SAVEPATH SAVEFILENAME],'w'); 
fprintf(fid,'%d\n',size(data,2)) ; % number of dimensions
fclose(fid);
dlmwrite([SAVEPATH SAVEFILENAME],data,'-append','delimiter',' '); % write rest of data

% run Klustakwik
temp2 = [KKPATH ' "' SAVEPATH Fname(temp(end)+1:end) '" 1'];
dos(temp2,'-echo');

CLUSFILENAME = [Fname(temp(end)+1:end) '.clu.1'];
a = load([SAVEPATH CLUSFILENAME]);
nk = a(1);
assign = a(2:end);



