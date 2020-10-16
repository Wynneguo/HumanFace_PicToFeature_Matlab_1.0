function GetLandmarks(Path_template,Path_template_lm,Dirpath_symmetry_txt,Dirpath_output)
%get lm of all
%input template.obj | template_lm_artifitial(27*4) | dir/samples
%output dir/samples_lm(27*3)
%导入template.obj template_lm %找不到对应点-需要找临近点再找索引
ptCloud_template=readOBJ(Path_template);
ptCloud_template=pointCloud(ptCloud_template);
template_lm=importdata(Path_template_lm);
% get lm
for i=1:size(SampleNames)
	LM_sample=zeros(27,3); 
	sample=importdata([Dirpath_symmetry_txt,SampleNames{i,1},'.txt']);
    ptCloud_sample=pointCloud(sample);
	for j=1:27
		point=template_lm(j,:);
		[indices,dists] = findNearestNeighbors(ptCloud_template,point(:,2:4),1); %找到17个landmark在点云中对应的点（距离为1的点）
		LM_sample(j,:)=ptCloud_sample.Location(indices,:);
		[m, ~] = size(LM_sample);
		fid=fopen([Dirpath_output,SampleNames{i,:},'.txt'], 'wt'); %按行写文件
		for k = 1 : m
			fprintf(fid, '%g\t', LM_sample(k,:));
			fprintf(fid, '\n');
        end
		fclose(fid);
    end
end
end