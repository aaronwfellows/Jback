function [hgoes,dgoes] = goesread(file)

hgoes=[]; %#ok<*NASGU>
dgoes=[];

if exist(file,'file')==2
  
   fid =  fopen(file,'r+','l');
   fseek(fid,0,'bof');

   newline=fgets(fid); 
   
else
  
  disp('no file')
  
end

icomma=findstr(newline,',');

a=[0 icomma length(newline)+1];
hgoes=[];

for j = 1:length(a)-1
    hgoes=strvcat(hgoes,newline(a(j)+1:a(j+1)-1));
end

% read data
dgoes = csvread(file,1,0)';

return


