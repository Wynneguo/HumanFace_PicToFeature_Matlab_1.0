function ChangeNoseFormat(SampleNames,Dirpath_Orinose,Dirpath_output)
% DirPath with '/' at the end
for i=1:size(SampleNames_1,1)
    nose=importdata([Dir_Orinose,num2str(SampleNames(i,1)),'.txt']);
    temp=split(nose(1,1),["(",",",")"]);
    nose=[temp(2,1),temp(3,1),temp(4,1)];
    temp=cellfun(@str2num,nose);
    nose=temp;
    [m,n]=size(nose);
    fid=fopen([Dir_output,num2str(SampleNames(i,1)),'.txt'],'wt');
    for j=1:n
        fprintf(fid,"%g",nose(1,j));
        fprintf(fid,' ');
    end
    fclose(fid);
end
end