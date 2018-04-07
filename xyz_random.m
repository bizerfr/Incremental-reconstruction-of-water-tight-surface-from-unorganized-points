clc;
clear;

[filename, pathname] = uigetfile('*.xyz', 'read .xyz file'); %Ñ¡ÔñÍ¼Æ¬ÎÄ¼þ
if isequal(filename,0)   
   msgbox('no selection');
else
   file_path=fullfile(pathname, filename);  %get the file path
end

d=importdata(file_path);

[rows,cols]=size(d);


id=randperm(rows);


fid=fopen([file_path(1:end-4),'_r.xyz'],'w');

t=zeros(1,3);
for i=1:rows
    t=d(id(i),:);

    
    fprintf(fid,'%f %f %f \r\n',t(1),t(2),t(3));
    
end

fclose(fid);