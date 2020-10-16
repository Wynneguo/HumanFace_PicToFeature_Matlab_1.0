function OrigObjMinusNosetip(SampleNames,Dirpath_origobj,Dirpath_nosetip,Dirpath_output)
for i=1:size(SampleNames,1)
    try
        [vertices,faces]=readOBJ([Dirpath_origobj,SampleNames{i,1},'_mesh.obj']);
    catch
        [vertices,faces]=readOBJ([Dirpath_origobj,SampleNames{i,1},'.obj']);
    end
    nosetip=load([Dirpath_nosetip,SampleNames{i,1},'.txt']);
    vertices_1=vertices-nosetip;
    vertface2obj(vertices_1,faces,[Dirpath_output,SampleNames{i,1},'.obj'])
end