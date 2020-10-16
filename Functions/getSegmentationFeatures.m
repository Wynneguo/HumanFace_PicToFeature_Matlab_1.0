function getSegmentationFeatures(SampleNames,Path_template,Path_wholefaceSImatrix,Dirpath_samples_txt,Mode_RVmatrix,n_levels,n_pc,Dirpath_output)
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
        saveas(gca,[Dirpath_output,'segmentation_pic\',num2str(k),'_',num2str(j),'.jpeg']);
    end
end

temp_label=LabelTotal(n_levels+1,:);
hf=scatter3(template_v(:,1),template_v(:,2),template_v(:,3),'filled','cdata',temp_label);
hf=gca; 
grid off;
axis off;
view([0 90]);
saveas(gca,[Dirpath_output,'segmentation_pic\ALL.jpeg']);

%% get PC
% save the indices of these several segments
%% extract PC
mkdir(Dirpath_output,'segmentation_indice\');
%LabelTotal=importdata("G:\2016_2018\segmentation_cov_label.mat");
a=0;
str1='%g\t';
for i=2:size(LabelTotal,1)
    for j=1:max(LabelTotal(i,:))
        indice=find(LabelTotal(i,:)==j);
        if ~isempty(indice)
            a=a+1; %total number of groups
            fid=fopen([Dirpath_output,'segmentation_indice\level_',num2str(i),'_',num2str(a),'.txt'],'wt');
            str1='%g\t';
            str2=repmat(str1,1,size(indice,1));
            fprintf(fid,str2,indice);
            fclose(fid);
        end
    end
end

SI_matrix=importdata(Path_wholefaceSImatrix);
[SI_matrix_s,MU_s,SIGMA_s]=zscore(SI_matrix);
Path_indices=[Dirpath_output,'segmentation_indice\'];
mkdir([Dirpath_output,'segmentation_pc\']);
dir2=[Dirpath_output,'segmentation_pc\'];
k=n_pc;
    
File=dir(fullfile(Path_indices,'*.txt'));
IndiceNames={File.name}';
IndiceNames=split(IndiceNames,'.');
IndiceNames=IndiceNames(:,1);

for i=1:size(IndiceNames,1)
    InNames=IndiceNames{i,:};
    indices=importdata([Path_indices,strcat(char(IndiceNames(i,:))),'.txt']);
    SI_matrix_i=SI_matrix_s(:,indices);
    %Nan
    %NaN process2
    SI_matrix_i(:,find(sum(isnan(SI_matrix_i))))=[];
    [coeff_SI,score_SI,latent_SI,~,explained_SI,mu_SI] = pca(SI_matrix_i);
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
        name(1,x+2)={['pc_',num2str(x)]};
    end
    A.Properties.VariableNames=name(1,:); 
    writetable(A,[dir2,'\',InNames,'_pc.txt'],'Delimiter',' ');
end
end



