function FormReflection(SampleNames,Dirpath_objminusnose,Dirpath_output)
% the reflection is the one before registration
for i=1:size(SampleNames,1)
    [obj.Vertices,obj.Faces]=readOBJ([Dirpath_objminusnose,SampleNames{i,1},'.obj']);
    obj.Vertices(:,1)=-obj.Vertices(:,1)
    vertface2obj(obj.Vertices,obj.Faces,[Dirpath_output,SampleNames{i,:},'.obj'])
end
end