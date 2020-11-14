function getSegmentationPics(SampleNames,Path_template,Dirpath_samples_txt,Mode_RVmatrix,n_levels,Dirpath_output)
%Mode_RVmatrix:cov/cor
%n_levels:num of levels of HierarchicalSegmentation(5)
[template_v,template_f]=readOBJ(Path_template);
ntps=size(template_v,1);
data=zeros(size(SampleNames,1),ntps*3);
tmp_line=zeros(1,ntps*3);

for i=1:size(SampleNames,1)
    vertex=importdata([Dirpath_samples_txt,SampleNames{i,:},'.txt']);
    [transform,vertex1]=procrustes(template_v,vertex);
    for j=1:size(vertex1,1)
        tmp_line(1,(1+3*(j-1)):3*j)=vertex1(j,:);
    end
    data(i,:)=tmp_line;   
end
%data=importdata("G:\2016_2018\data_for_RV.txt"); %矩阵过大无法打开

RV_matrix=getRVmatrix(data,Mode_RVmatrix);
LabelTotal=HierarchicalSegmentationUsingSpectralClustering(RV_matrix,n_levels);
save([Dirpath_output,'segmentation_',Mode_RVmatrix,'_label.mat'],'LabelTotal');
%LabelTotal=importdata([Dirpath_output,'segmentation_',Mode_RVmatrix,'_label.mat']);

%% Draw segments
mkdir([Dirpath_output,'segmentation_pic\']);
for k=1:size(LabelTotal,1)
    temp_label=LabelTotal(k,:);
    uni=unique(temp_label);
    for j=uni
        x=zeros(size(temp_label));
        i=find(temp_label==j);
        x(i)=temp_label(i);
        temp_label_1=x;
        colormap('Summer');
        hf=scatter3(template_v(:,1),template_v(:,2),template_v(:,3),'filled','cdata',temp_label_1);
        hf=gca; 
        grid off;
        axis off;
        view([0 90]);
        saveas(gca,[Dirpath_output,'segmentation_pic\',num2str(k),'_',num2str(j),'_level_',num2str(n_levels),'.jpeg']);
    end
end

temp_label=LabelTotal(n_levels+1,:);
hf=scatter3(template_v(:,1),template_v(:,2),template_v(:,3),'filled','cdata',temp_label);
hf=gca; 
grid off;
axis off;
view([0 90]);
saveas(gca,[Dirpath_output,'segmentation_pic\','level_',num2str(n_levels),'_ALL.jpeg']);
disp(['finish level：' num2str(k)]);
end
