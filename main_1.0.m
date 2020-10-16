
%% registration
%% input id | original_obj(obj_1) | original_nosetip(nosetip_1)

fid=fopen("D:\guolu\id\6004id_pic.txt");
SampleNames=textscan(fid,'%s');
fclose(fid);
SampleNames=SampleNames{1,1};

%% change the fomat of nosetips
SampleNames_1=importdata("C:\Users\dell\Desktop\613samples_20200210_imputation_filterINFO.txt");
SampleNames_1=SampleNames_1(:,1);
for i=1:size(SampleNames_1,1)
    nose=importdata(['D:\guolu\nosetip_1\',num2str(SampleNames_1(i,1)),'.txt']);
    temp=split(nose(1,1),["(",",",")"]);
    nose=[temp(2,1),temp(3,1),temp(4,1)];
    temp=cellfun(@str2num,nose);
    nose=temp;
    [m,n]=size(nose);
    fid=fopen(['D:\guolu\temp\',num2str(SampleNames_1(i,1)),'.txt'],'wt');
    for j=1:n
        fprintf(fid,"%g",nose(1,j));
        fprintf(fid,' ');
    end
    fclose(fid);
end
    
%% show the nosetip || make sure it is the right match between image and nosetip
[vertices,faces]=readOBJ('D:\guolu\obj_1\0005.obj');
pcshow(vertices);
hold on;
nosetip=load('D:\guolu\nosetip_1\0005.txt');
scatter3(nosetip(1,1),nosetip(1,2),nosetip(1,3));

%% orig_obj minus nosetip || convert the cordination of nosetip into (0,0)
for i=1:size(SampleNames,1)
    try
        [vertices,faces]=readOBJ(['D:\guolu\obj_1\',SampleNames{i,1},'_mesh.obj']);
    catch
        [vertices,faces]=readOBJ(['D:\guolu\obj_1\',SampleNames{i,1},'.obj']);
    end
    nosetip=load(['D:\guolu\nosetip_1\',SampleNames{i,1},'.txt']);
    vertices_1=vertices-nosetip;
    vertface2obj(vertices_1,faces,['D:\guolu\obj_minusnose\',SampleNames{i,1},'.obj'])
end

%% Formation of reflection || symmetry=(reflection+individual)/2
SampleNames=regis_object;
Path='D:\guolu\obj_minusnose\';
for i=1:size(SampleNames,1)
    [obj.Vertices,obj.Faces]=readOBJ([Path,SampleNames{i,1},'.obj']);
    %Vertices=array2table(obj.Vertices);
    %writetable(Vertices,['D:\guolu\txt_individual\',SampleNames{i,1},'.txt'],'Delimiter',' ');
    obj.Vertices(:,1)=-obj.Vertices(:,1)
    vertface2obj(obj.Vertices,obj.Faces,['D:\guolu\temp_reflection615\',SampleNames{i,:},'.obj'])
end

%% rigid_registration || just an alignment with one chosen pic
%input dir/samples_obj | template_DefinedFromSamples_obj | dir/samples_lm_true/ (TemplateIncluded)
%output dir/samples_afterRigidRegis_obj/

addpath(genpath('D:\guolu\meshmonk')) ;
studypath = 'D:\guolu\obj_minusnose';
cd(studypath);

for i=1:size(SampleNames,1)
	%%load template
	Template=shape3D;
	filename='template.obj';     
	path='D:/guolu/';
	importWavefront(Template,filename,path,[]) %%导入obj func-importWavefront
    demoFace = shape3D;  %==>demoFace.Type=shape 3D
    try
        filename =[SampleNames{i,:},'.obj']; 
    catch 
        filename =[SampleNames{i,:},'_mesh.obj'];
    end
    path='D:\guolu\obj_minusnose\';
    importWavefront(demoFace,filename,path,[]);
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
    vertface2obj(obj.FloatingShape.Vertices,obj.FloatingShape.Faces,['D:\guolu\obj_individual\',SampleNames{i,:},'.obj'])
    %x=array2table(obj.FloatingShape.Vertices);
    %writetable(x,['D:\guolu\txt_reflection\',SampleNames{i,:},'.txt'],'Delimiter',' ','WriteVariableNames',false);

end
%% Formation of Symmetry
for i=1:size(SampleNames,1)
    [v1,f]=readOBJ(['D:\guolu\aaaa_diffobj\diffobj_individual\',SampleNames{i,:},'.obj']);
    fid=fopen(['D:\guolu\aaaa_diffobj\txt_reflection\',SampleNames{i,:},'.txt']);
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
    vertface2obj(v,f,['D:\guolu\aaaa_diffobj\symmetryOBJ\',SampleNames{i,:},'.obj'])
    v=array2table(v);
    writetable(v,['D:\guolu\aaaa_diffobj\symmetryTXT\',SampleNames{i,:},'.txt'],'Delimiter',' ','WriteVariableNames',false);
end
%% get lm of all
%input template.obj | template_lm_artifitial(fromStep1) | dir/samples
%output dir/samples_lm
%导入template.obj template_lm %找不到对应点-需要找临近点再找索引
ptCloud_template=readOBJ('D:\guolu\template.obj');
ptCloud_template=pointCloud(ptCloud_template);
template_lm=importdata('D:\guolu\template_lm.txt');
% get lm
for i=1:size(SampleNames)
	LM_sample=zeros(27,3); 
	sample=importdata(['D:\guolu\symmetryTXT\',SampleNames{i,1},'.txt']);
    ptCloud_sample=pointCloud(sample);
	for j=1:27
		point=template_lm(j,:);
		[indices,dists] = findNearestNeighbors(ptCloud_template,point(:,2:4),1); %找到17个landmark在点云中对应的点（距离为1的点）
		LM_sample(j,:)=ptCloud_sample.Location(indices,:);
		[m, ~] = size(LM_sample);
		fid=fopen(['D:\guolu\landmarks\',SampleNames{i,:},'.txt'], 'wt'); %按行写文件
		for k = 1 : m
			fprintf(fid, '%g\t', LM_sample(k,:));
			fprintf(fid, '\n');
        end
		fclose(fid);
    end
end

%% The Get features Now!
Path_of_Samplelm='D:\guolu\landmarks\';
Path_of_features='D:\guolu\features\';
template_lm=importdata('D:\guolu\template_lm.txt');
lm_mark=template_lm(:,1);%the mark of points
num_of_points=size(template_lm,1);
PointsCombinations=combntns([1:num_of_points],2); % combination of all the distance conditions
sample_dist={};

%get features
for k=1:size(PointsCombinations,1)
    sample_lm_matirx=zeros(2,3);
    cp_1=PointsCombinations(k,1);
    cp_2=PointsCombinations(k,2);
    for i=1:size(SampleNames,1)
        sample_lm=importdata([Path_of_Samplelm,strcat(char(SampleNames(i,:))),'.txt']);
        [~,sample_lm]=procrustes(template_lm(:,2:4),sample_lm); %Procrustes the lm(translation, reflection, orthogonal rotation, and scaling)
        sample_lm=[lm_mark,sample_lm];
        for j=1:size(sample_lm,1)
            if sample_lm(j,1)==cp_1
                sample_lm_matirx(1,:)=sample_lm(j,2:4);
            end
            if sample_lm(j,1)==cp_2
                sample_lm_matirx(2,:)=sample_lm(j,2:4);
            end
        end
        dis=pdist(sample_lm_matirx);
        if k==1
            sample_dist{i,1}=SampleNames{i,:};
        end
        sample_dist{i,1+k}=dis;
    end
end
%change the colnames
sample_dist=cell2table(sample_dist); 
sample_dist.Properties.VariableNames{1,1}='ID';
for y=1:size(PointsCombinations,1)
    cp_1=PointsCombinations(y,1);
    cp_2=PointsCombinations(y,2);      
    sample_dist.Properties.VariableNames{1,y+1}=[num2str(cp_1),'_',num2str(cp_2)];  
end
%save distance features
writetable(sample_dist,[Path_of_features,'features_distance.txt'],'Delimiter',' ','WriteVariableNames',true);
%a=readtable([Path_of_features,'features_distance.txt'],'Format',['%s',repmat('%f',1,351)]);
%% get angle features
AngleFeatures={}
%Rf1 眼睛大小/眼间距
Rf1={};
lm_withindex=[];
for i=1:size(SampleNames,1)
    lm=importdata([Path_of_Samplelm,SampleNames{i,:},'.txt']);
    lm_withindex=[lm_mark,lm];
    for j=1:size(lm_withindex,1)
        if lm_withindex(j,1)==7
            a=lm_withindex(j,2:4);
        elseif lm_withindex(j,1)==3
            b=lm_withindex(j,2:4);
        elseif lm_withindex(j,1)==4
            c=lm_withindex(j,2:4);
        elseif lm_withindex(j,1)==9
            d=lm_withindex(j,2:4);
        end
    end
    ratio=((norm(a-b)+norm(c-d))/2)/norm(b-c);
    Rf1{i,1}=SampleNames{i,:};
    Rf1{i,2}=ratio;
end
AngleFeatures=cat(2,AngleFeatures,Rf1);

% Rf2 上眼眶-鼻底/鼻底-下巴
Rf2={};
lm_withindex=[];
for i=1:size(SampleNames,1)
    lm=importdata([Path_of_Samplelm,SampleNames{i,:},'.txt']);
    lm_withindex=[lm_mark,lm];
    for j=1:size(lm_withindex,1)
        if lm_withindex(j,1)==15
            a=lm_withindex(j,3);
        elseif lm_withindex(j,1)==10
            b=lm_withindex(j,3);
        elseif lm_withindex(j,1)==16
            c=lm_withindex(j,3);
        elseif lm_withindex(j,1)==11
            d=lm_withindex(j,3);
        end
    end
    ratio=((norm(a-b)+norm(c-b))/2)/norm(b-d);
    Rf2{i,1}=strcat(char(SampleNames(i,:)));
    Rf2{i,2}=ratio;
end
AngleFeatures=cat(2,AngleFeatures,Rf2(:,2));
% Anglef1
Af1={};
lm_withindex=[];
for i=1:size(SampleNames,1)
    lm=importdata([Path_of_Samplelm,SampleNames{i,:},'.txt']);
    lm_withindex=[lm_mark,lm];
    for j=1:size(lm_withindex,1)
        if lm_withindex(j,1)==17
            a=lm_withindex(j,2:4);
        elseif lm_withindex(j,1)==18
            b=lm_withindex(j,2:4);
        elseif lm_withindex(j,1)==10
            c=lm_withindex(j,2:4);
        end
    end
    l1=b-a;
    l2=b-c;
    angle=acos(l1*l2'/(norm(l1)*norm(l2)));
    Af1{i,1}=strcat(char(SampleNames(i,:)));
    Af1{i,2}=angle;
end
AngleFeatures=cat(2,AngleFeatures,Af1(:,2));
% Anglef2
Af2={};
lm_withindex=[];
for i=1:size(SampleNames,1)
    lm=importdata([Path_of_Samplelm,SampleNames{i,:},'.txt']);
    lm_withindex=[lm_mark,lm];
    for j=1:size(lm_withindex,1)
        if lm_withindex(j,1)==18
            a=lm_withindex(j,2:4);
        elseif lm_withindex(j,1)==10
            b=lm_withindex(j,2:4);
        elseif lm_withindex(j,1)==21
            c=lm_withindex(j,2:4);
        end
    end
    l1=b-a;
    l2=b-c;
    angle=acos(l1*l2'/(norm(l1)*norm(l2)));
    Af2{i,1}=strcat(char(SampleNames(i,:)));
    Af2{i,2}=angle;
end
AngleFeatures=cat(2,AngleFeatures,Af2(:,2));

% Af3
Af3={};
lm_withindex=[];
for i=1:size(SampleNames,1)
    lm=importdata([Path_of_Samplelm,SampleNames{i,:},'.txt']);
    lm_withindex=[lm_mark,lm];
    for j=1:size(lm_withindex,1)
        if lm_withindex(j,1)==10
            a=lm_withindex(j,2:4);
        elseif lm_withindex(j,1)==21
            b=lm_withindex(j,2:4);
        elseif lm_withindex(j,1)==2
            c=lm_withindex(j,2:4);
        end
    end
    l1=b-a;
    l2=b-c;
    angle=acos(l1*l2'/(norm(l1)*norm(l2)));
    Af3{i,1}=strcat(char(SampleNames(i,:)));
    Af3{i,2}=angle;
end
AngleFeatures=cat(2,AngleFeatures,Af3(:,2));
%  Af4
Af4={};
lm_withindex=[];
for i=1:size(SampleNames,1)
    lm=importdata([Path_of_Samplelm,SampleNames{i,:},'.txt']);
    lm_withindex=[lm_mark,lm];
    for j=1:size(lm_withindex,1)
        if lm_withindex(j,1)==21
            a=lm_withindex(j,2:4);
        elseif lm_withindex(j,1)==2
            b=lm_withindex(j,2:4);
        elseif lm_withindex(j,1)==11
            c=lm_withindex(j,2:4);
        end
    end
    l1=b-a;
    l2=b-c;
    angle=acos(l1*l2'/(norm(l1)*norm(l2)));
    Af4{i,1}=strcat(char(SampleNames(i,:)));
    Af4{i,2}=angle;
end
AngleFeatures=cat(2,AngleFeatures,Af4(:,2));
% Af5
Af5={};
lm_withindex=[];
for i=1:size(SampleNames,1)
    lm=importdata([Path_of_Samplelm,SampleNames{i,:},'.txt']);
    lm_withindex=[lm_mark,lm];
    for j=1:size(lm_withindex,1)
        if lm_withindex(j,1)==19
            a=lm_withindex(j,2:4);
        elseif lm_withindex(j,1)==17
            b=lm_withindex(j,2:4);
        elseif lm_withindex(j,1)==20
            c=lm_withindex(j,2:4);
        end
    end
    l1=b-a;
    l2=b-c;
    angle=acos(l1*l2'/(norm(l1)*norm(l2)));
    Af5{i,1}=strcat(char(SampleNames(i,:)));
    Af5{i,2}=angle;
end
AngleFeatures=cat(2,AngleFeatures,Af5(:,2));
% Af6
Af6={};
lm_withindex=[];
for i=1:size(SampleNames,1)
    lm=importdata([Path_of_Samplelm,SampleNames{i,:},'.txt']);
    lm_withindex=[lm_mark,lm];
    for j=1:size(lm_withindex,1)
        if lm_withindex(j,1)==19
            a=lm_withindex(j,2:4);
        elseif lm_withindex(j,1)==18
            b=lm_withindex(j,2:4);
        elseif lm_withindex(j,1)==20
            c=lm_withindex(j,2:4);
        end
    end
    l1=b-a;
    l2=b-c;
    angle=acos(l1*l2'/(norm(l1)*norm(l2)));
    Af6{i,1}=strcat(char(SampleNames(i,:)));
    Af6{i,2}=angle;
end
AngleFeatures=cat(2,AngleFeatures,Af6(:,2));
% Area nose
Area={};
lm_withindex=[];
for i=1:size(SampleNames,1)
    lm=importdata([Path_of_Samplelm,SampleNames{i,:},'.txt']);
    lm_withindex=[lm_mark,lm];
    for j=1:size(lm_withindex,1)
        if lm_withindex(j,1)==17
            a=lm_withindex(j,2:4);
        elseif lm_withindex(j,1)==19
            b=lm_withindex(j,2:4);
        elseif lm_withindex(j,1)==10
            c=lm_withindex(j,2:4);
        elseif lm_withindex(j,1)==20
            d=lm_withindex(j,2:4);
        end
    end
    l1=norm(a-b);
    l2=norm(b-c);
    l3=norm(c-a);
    
    l4=norm(d-a);
    l5=norm(c-d);
    p1=(l1+l2+l3)/2 ;
    S1=sqrt(p1*(p1-l1)*(p1-l2)*(p1-l3));
    p2=(l4+l5+l3)/2 ;
    S2=sqrt(p2*(p2-l4)*(p2-l5)*(p2-l3));
    S=S1+S2;
    Area{i,1}=strcat(char(SampleNames(i,:)));
    Area{i,2}=S;
end
AngleFeatures=cat(2,AngleFeatures,Area(:,2));
%change the colnames
AngleFeatures=cell2table(AngleFeatures); 
name={'Rf1','Rf2','Af1','Af2','Af3','Af4','Af5','Af6','Area'};
AngleFeatures.Properties.VariableNames{1,1}='ID';
for i=1:size(name,2)
    AngleFeatures.Properties.VariableNames{1,i+1}=name{1,i};
end
%save angle features
writetable(AngleFeatures,[Path_of_features,'AngleFeatures.txt'],'Delimiter',' ','WriteVariableNames',true);
%a=readtable([Path_of_features,'features_distance.txt'],'Format',['%s',repmat('%f',1,351)]);

%% get SI (Whole)
%得到所有samples指定数量的pc [FID IID PC1 PC2 PCn] (SI Accessible|GS invalid)
%input template.ply | Landmark(帮助得到合适的脸部范围_轮廓) |sample.id.txt 
%input dir/sample.obj/  (after_NRigid_individual.obj)
%output integration_face_SI

% import data
ptCloud_template=readOBJ('D:\guolu\template.obj');
ptCloud_template=pointCloud(ptCloud_template);
LM_data=importdata('D:\guolu\template_lm.txt');
LM_template=LM_data(1:27,2:4);
% draw district
% % coordinate=ptCloud_template.Location(indices,:);
% % pcshow(coordinate);

% get curvature and shape_index
ntps=size(ptCloud_template.Location,1);
SI_matrix=zeros(size(SampleNames,1),ntps);
obj_Path='D:\guolu\symmetryOBJ\';
for i=1:size(SampleNames,1)
    options.curvature_smoothing = 10;
    options.verb = 0;
    SI=[];
    Cgauss=[];
    [vertex,faces]=readOBJ([obj_Path,strcat(char(SampleNames(i,:))),'.obj']);
    [Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(vertex,faces,options); 
    for j = 1:ntps
       SI(j,1)=(-2/pi)*atan((Cmax(j,1)+Cmin(j,1))/(Cmax(j,1)-Cmin(j,1)));   %shape index
    end
    SI_face=SI; 
    SI_matrix(i,:)=SI_face';
    disp(i);
end
%%
%Storage of the matrix
savepath='D:\guolu\features\'
[m, n] = size(SI_matrix);
fid=fopen([savepath,'SI_matrix_intergral_face_original.txt'], 'wt');
for i = 1 : m
	fprintf(fid, '%g\t', SI_matrix(i,:));
	fprintf(fid, '\n');
end
fclose(fid);
%% NaN casuse the error in PCA
%NaN process1(replace Nan with colmean with at least one real num in this
%col)
for j=1:size(SI_matrix,1)
    for k=1:size(SI_matrix,2)
        if isnan(SI_matrix(j,k))
            SI_matrix(j,k)=nanmean(SI_matrix(:,k)); 
        end
    end
end
%NaN process2
SI_matrix_1=SI_matrix;
SI_matrix_1(:,find(sum(isnan(SI_matrix))))=[];
%NaN process3(not be used for the time being)
col=find(sum(isnan(SI_matrix_1)));
col(2,:)=0;
col_orig=1:size(vertex,1);
col_diff=setdiff(col_orig,col);
for i=1:size(col,2)
    while col(2,i)~=0
        j=1;
        bl=[sum(col_diff==col(1,i)-j)==1 ,sum(col_diff==col(1,i)+j)==1];
        if bl==[1,0]
            col(2,i)=col(1,i)-j;
        elseif bl==[0,1]
            col(2,i)=col(1,i)+j;
        else
            j=j+1;
        end
    end
end

%% 
%read saved matrix data
SI_matrix=importdata([savepath,'SI_matrix_intergral_face.txt']);
% Standardization **
[SI_matrix_s,MU_s,SIGMA_s]=zscore(SI_matrix_1);
% get scores & PC
[coeff_SI,score_SI,latent_SI,~,explained_SI,mu_SI] = pca(SI_matrix_s);
% merge pheno file
k=50;
A={};         
A(:,1)=SampleNames(:,1);  
A(:,2)=SampleNames(:,1);  
for x=1:k
    for i=1:size(SampleNames,1)
        A(i,x+2)={num2str(score_SI(i,x))};
    end
end
A=cell2table(A); 
name={};
name(1,1)={'FID'};
name(1,2)={'IID'};
for x=1:k
    name(1,x+2)={['face_SI_',num2str(x)]};
end
for i=1:size(name,2)
    A.Properties.VariableNames{1,i}=name{1,i}; 
end
writetable(A,[savepath,'integration_face_SI.txt'],'Delimiter',' ');

%% get SI(nose)
% % %% get the exact region(skip_replate by nose template)
% % % import data
% % ptCloud_template= template.Vertices;
% % Ntp1 = LM_template(10,:);  %nose_lm 
% % point_c=(LM_template(2,:)+LM_template(9,:))/2;%计算半径
% % radius=ceil(pdist([LM_template(2,:);point_c]))+15%区域较大一些
% % indices = findNeighborsInRadius(ptCloud_template,Ntp1,radius);     %find points within a region of interest.
%% get SI of nose
PartName='Nose';
Path_template='D:/guolu/template.obj';
Path_template_part='D:/guolu/template_nose_from geomagic.obj';
Dirpath_symmetryobj='D:\guolu\symmetryOBJ\';
Dirpath_ouput='D:\guolu\features\';
n_pc=50;
GetSIfeatures(SampleNames,PartName,Path_template,Path_template_part,Dirpath_symmetryobj,Dirpath_ouput,n_pc)
%% get SI of segmentation
Path_wholefaceSImatrix='D:\guolu\features\SI_matrix_intergral_face_original.txt';
Dirpath_samples_txt='D:\guolu\symmetryTXT\';
Mode_RVmatrix='cov';
n_levels=5;
n_pc=50;
Dirpath_output='D:\guolu\features\';
getSegmentationFeatures(SampleNames,Path_template,Path_wholefaceSImatrix,Dirpath_samples_txt,Mode_RVmatrix,n_levels,n_pc,Dirpath_output)




