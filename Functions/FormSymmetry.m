function FormSymmetry(SampleNames,Dirpath_individual,Dirpath_reflection_txt,Dirpath_output,Dirpath_ouput_txt)
   %Dirpath_individual&Dirpath_reflection involved Afterregistered pic
for i=1:size(SampleNames,1)
    [v1,f]=readOBJ([Dirpath_individual,SampleNames{i,:},'.obj']);
    fid=fopen([Dirpath_reflection,SampleNames{i,:},'.txt']);
    v2=textscan(fid,'%s');
    fclose(fid);
    v2=v2{1,1};
    for j=1:size(v2,1)/3
        temp(j,1)=str2num(v2{3*(j-1)+1,1});
        temp(j,2)=str2num(v2{3*(j-1)+2,1});
        temp(j,3)=str2num(v2{3*(j-1)+3,1});
    end
    v2=temp;
    v=(v1+v2)/2;
    vertface2obj(v,f,[Dirpath_output,SampleNames{i,:},'.obj'])
    v=array2table(v);
    writetable(v,[Dirpath_ouput_txt,SampleNames{i,:},'.txt'],'Delimiter',' ','WriteVariableNames',false);
end
end
