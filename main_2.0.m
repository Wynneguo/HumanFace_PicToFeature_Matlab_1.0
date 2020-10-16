%% 20201006 Pic Process 
%% Input ID 
fid=fopen("D:\guolu\id\6004id_pic.txt");
SampleNames=textscan(fid,'%s');
fclose(fid);
SampleNames=SampleNames{1,1};
%%  change the fomat of nosetips
Dirpath_orinose='';
Dirpath_nose='';
ChangeNoseFormat[SampleNames,Dirpath_orinose,Dirpath_nose];
%% show the nosetip || make sure it is the right match between image and nosetip
[vertices,faces]=readOBJ('');
pcshow(vertices);
hold on;
nosetip=load('');
scatter3(nosetip(1,1),nosetip(1,2),nosetip(1,3));
%% orig_obj minus nosetip || convert the cordination of nosetip into (0,0)
Dirpath_origobj='';
Dirpath_nosetip='';
Dirpath_oriminusnose='';
OrigObjMinusNosetip[SampleNames,Dirpath_origobj,Dirpath_nosetip,Dirpath_oriminusnose];
%% Formation of reflection || symmetry=(reflection+individual)/2
Dirpath_MinusNoseReflection='';
FormReflection[SampleNames,Dirpath_objminusnose,Dirpath_MinusNoseReflection];
%% Registration 
addpath(genpath('D:\guolu\meshmonk')) ;
TemplatePath='';
Dirpath_oriminusnose='';
Dirpath_MinusNoseReflection='';
Dirpath_individual='';
Dirpath_reflection='';
Dirpath_txt_individual='';
Dirpath_txt_reflection='';
Registration_comflation[SampleNames,TemplatePath,Dirpath_oriminusnose,Dirpath_individual,Dirpath_txt_individual];
Registration_comflation[SampleNames,TemplatePath,Dirpath_MinusNoseReflection,Dirpath_reflection,Dirpath_txt_reflection];
%% Formation of Symmetry
Dirpath_individual='';
Dirpath_reflection_txt='';
Dirpath_output='';
Dirpath_ouput_txt='';
FormSymmetry[SampleNames,Dirpath_individual,Dirpath_reflection_txt,Dirpath_output,Dirpath_ouput_txt]
%% Get landmarks of all
Path_template='/template.obj';
Path_template_lm='/template_lm.txt';
Dirpath_symmetry_txt='';
Dirpath_lm='';
GetLandmarks[Path_template,Path_template_lm,Dirpath_symmetry_txt,Dirpath_lm]
%% Get features
Path_template_lm='';
Dirpath_of_Samplelm='';
Dirpath_features='';
Getdistfeatures[SampleNames,Path_template_lm,Dirpath_of_Samplelm,Dirpath_features];
GetAngleFeatures[SampleNames,Path_template_lm,Dirpath_of_Samplelm,Dirpath_features];
%a=readtable([Path_of_features,'features_distance.txt'],'Format',['%s',repmat('%f',1,351)]);

PartName='Nose';
Path_template='D:/guolu/template.obj'
Path_template_part='D:/guolu/template_nose_from geomagic.obj';
Dirpath_symmetryobj='D:\guolu\symmetryOBJ\';
Dirpath_ouput='D:\guolu\features\';
n_pc=50;
GetSIfeatures(SampleNames,PartName,Path_template,Path_template_part,Dirpath_symmetryobj,Dirpath_ouput,n_pc)

Path_wholefaceSImatrix='D:\guolu\features\SI_matrix_intergral_face_original.txt';
Dirpath_samples_txt='D:\guolu\symmetryTXT\';
Mode_RVmatrix='cov';
n_levels=5;
n_pc=50;
Dirpath_output='D:\guolu\features\';
getSegmentationFeatures(SampleNames,Path_template,Path_wholefaceSImatrix,Dirpath_samples_txt,Mode_RVmatrix,n_levels,n_pc,Dirpath_output)



