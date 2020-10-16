function Registration_comflation(SampleNames,TemplatePath,Dirpath_inputobj,Dirpath_output,Dirpath_output_txt)
%combine the two registration process (rigid&nonrigid)
%input dir/samples_obj | template_DefinedFromSamples_obj 
%output dir/samples_afterRigidRegis_obj/

for i=1:size(SampleNames,1)
	%%load template
	Template=shape3D;
	filename='template.obj';     
	importWavefront(Template,filename,TemplatePath,[]) %%导入obj func-importWavefront
    demoFace = shape3D;  %==>demoFace.Type=shape 3D
    try
        filename =[SampleNames{i,:},'.obj']; 
    catch 
        filename =[SampleNames{i,:},'_mesh.obj'];
    end
    importWavefront(demoFace,filename,Dirpath_inputobj,[]);
    % demoFace.FlipNormals=1; %Only for reflection
    % Mapping demo face
    % Create an instance of the Shape Mapper class (obj)
    obj = ShapeMapper;
    obj.FloatingShape = clone(Template);      % Assign the Template to floating shape
    obj.TargetShape = clone(demoFace);        % Assign the face to map to target shape
    %% Rigid Registration
    obj.Display = false;                       % If false does not displays mapping process. When true, two windows appear:
                                              % One shows the alignment of the template (white mesh) onto the target. The other shows how the
                                              % template is being transformed during the registration, with yellow parts of the mesh
                                              % classified as inliers and blue parts of the mesh classified as outliers. 
                                              %两幅图 分别展示了template向target alignment-白网&&template是如何变形的  黄色部分 inliers/蓝色部分 outliers

    % Predefined rigid mapping parameters, can be adjusted as desired mapping的各种参数-可调 
    obj.NumIterations = 80;                   % Number of iterations for the non-rigid registration 迭代次数
    obj.InlierKappa = 3;                      % Threshold to consider a point as an outlier (+/- k times the standard deviation) 离群点阈值
    obj.TransformationType = 'rigid';         
    obj.UseScaling = true;                    
    obj.CorrespondencesNumNeighbours = 3;     % Number of k-nearest neighbors  k近邻个数
    obj.CorrespondencesFlagThreshold = 0.9;
    obj.CorrespondencesSymmetric = true;      % If false：push forces are calculated (typical one-to-one correspondences). 
                                              % When true: Pull-and-push forces are calculated ???
    obj.CorrespondencesEqualizePushPull = false;
    obj.InlierUseOrientation = true;
    obj.FlagFloatingBoundary = true;
    obj.FlagTargetBoundary = true;
    obj.FlagTargetBadlySizedTriangles = true; % Pruning large sized triangles  对大的三角形区域的删除修剪？下为‘大’的定义
    obj.TriangleSizeZscore = 6;               % A large triangle is considered to fall x times the standard deviation from the mean
    obj.UpSampleTarget = false;
    % Mapping
    tic;map(obj);time = toc;  %执行函数时间记录 tic Operation toc
    fprintf ( SampleNames{i,:}, '  The Rigid Registration took %f seconds to run.\n', time );
 
    %% Non-Rigid Registration
    obj.Display = true;                       % Same as above
    % Predefined rigid mapping parameters, can be adjusted as desired 
    nr = 200;                                 % Number of iterations for the non-rigid registration
    obj.NumIterations = nr;                   % Number of iterations for the non-rigid registration
    obj.TransformNumNeighbors = 80;           % Number of neighbors regularization
    obj.CorrespondencesNumNeighbours = 3;     % Number of k-nearest neighbors
    obj.CorrespondencesSymmetric = true;      % If false push forces are calculated (typical one-to-one correspondences). 
                                              % When true: Pull-and-push forces are calculated 
    obj.CorrespondencesFlagThreshold = 0.9;
    obj.CorrespondencesUseOrientation = true;
    obj.CorrespondencesEqualizePushPull =false;
    obj.InlierKappa = 12;                     % Threshold to consider a point as an outlier (+/- k times the standard deviation)
    obj.InlierUseOrientation = true;
    obj.FlagFloatingBoundary = true;
    obj.FlagTargetBoundary = true;
    obj.FlagTargetBadlySizedTriangles = true; % Pruning large sized triangles 
    obj.TriangleSizeZscore = 6;               % A large triangle is considered to fall x times the standard deviation from the mean
    obj.UpSampleTarget = false;
    obj.UseScaling = 1;
    obj.TransformationType = 'nonrigid';
    obj.TransformSigma = 3;
    obj.TransformNumViscousIterationsStart = nr;
    obj.TransformNumViscousIterationsEnd = 1;
    obj.TransformNumElasticIterationsStart = nr;
    obj.TransformNumElasticIterationsEnd = 1;
    % Mapping
    tic;map(obj);time = toc;
    fprintf ( SampleNames{i,:}, '  The Non-Rigid Registration took %f seconds to run.\n', time );
%     % Visualize Mapping result
%     vref = viewer3D;                         % Create an instace of the viewer
%     obj.TargetShape.ViewMode = 'points';     % Display as a point cloud
%     obj.FloatingShape.ViewMode = 'points';
%     obj.TargetShape.VertexSize = 10;         % Setup the size of the points
%     obj.FloatingShape.VertexSize = 12;
%     viewer(obj.FloatingShape,vref);          % Display the  Shape in the viewer
%     viewer(obj.TargetShape,vref);
    vertface2obj(obj.FloatingShape.Vertices,obj.FloatingShape.Faces,[Dirpath_output,SampleNames{i,:},'.obj'])
    x=array2table(obj.FloatingShape.Vertices);
    writetable(x,[Dirpath_output_txt,SampleNames{i,:},'.txt'],'Delimiter',' ','WriteVariableNames',false);

end
end