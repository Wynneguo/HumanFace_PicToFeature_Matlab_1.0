function GetSIfeatures_part(SampleNames,PartName,Path_SI_matrix_intergral_face,Path_template,Path_template_part,Dirpath_ouput,n_pc)
%get the exact pc [FID IID PC1 PC2 PCn] (SI Accessible|GS invalid)
%premise:Get the original whole face SI_matrix first
%partname:Nose/Face/Others 
%Landmark(for limiting the spheres of face_Did no use temporarily) 
%output SI.txt & PC.txt
SI_matrix=importdata(Path_SI_matrix_intergral_face);
% import data
ptCloud_template=readOBJ(Path_template);
ptCloud_template=pointCloud(ptCloud_template);
ptCloud_template_part=readOBJ(Path_template_part);
ptCloud_template_part=pointCloud(ptCloud_template_part);

%LM_data=importdata('D:\guolu\template_lm.txt');
%LM_template=LM_data(1:27,2:4);
% draw district
% % coordinate=ptCloud_template.Location(indices,:);
% % pcshow(coordinate);

% get curvature and shape_index
part_indices=[];
for i=1:size(ptCloud_template_part.Location,1)
    [indices,dists] = findNearestNeighbors(ptCloud_template,ptCloud_template_part.Location(i,:),1);
    part_indices=[part_indices;indices];
end
ntps=size(part_indices,1);
SI_matrix=SI_matrix(:,part_indices);

%%
%Storage of the matrix
savepath=Dirpath_ouput
[m, n] = size(SI_matrix);
fid=fopen([savepath,'SI_matrix_',PartName,'.txt'], 'wt');
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
% % %NaN process3(not be used for the time being)
% % col=find(sum(isnan(SI_matrix_1)));
% % col(2,:)=0;
% % col_orig=1:size(vertex,1);
% % col_diff=setdiff(col_orig,col);
% % for i=1:size(col,2)
% %     while col(2,i)~=0
% %         j=1;
% %         bl=[sum(col_diff==col(1,i)-j)==1 ,sum(col_diff==col(1,i)+j)==1];
% %         if bl==[1,0]
% %             col(2,i)=col(1,i)-j;
% %         elseif bl==[0,1]
% %             col(2,i)=col(1,i)+j;
% %         else
% %             j=j+1;
% %         end
% %     end
% % end

%% 
%read saved matrix data
%SI_matrix=importdata([savepath,'.txt']);
% Standardization **
[SI_matrix_s,MU_s,SIGMA_s]=zscore(SI_matrix_1);
% get scores & PC
[coeff_SI,score_SI,latent_SI,~,explained_SI,mu_SI] = pca(SI_matrix_s);
% merge pheno file
k=n_pc;
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
    name(1,x+2)={[PartName,'_SI_',num2str(x)]};
end
for i=1:size(name,2)
    A.Properties.VariableNames{1,i}=name{1,i}; 
end
writetable(A,[savepath,PartName,'_SI_PC_',n_pc,'.txt'],'Delimiter',' ');
end
