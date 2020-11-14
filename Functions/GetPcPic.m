%For part which have ptcloud_template
function GetPcPic(Path_template_part,Path_part_SImatrix,n_pc,Dirpath_output,PartName)
ptCloud_template=readOBJ(Path_template_part);
SI_matrix=importdata(Path_part_SImatrix);
[SI_matrix_i,MU_s,SIGMA_s]=zscore(SI_matrix);

Nan_col=find(sum(isnan(SI_matrix_i)));
SI_matrix_i(:,find(sum(isnan(SI_matrix_i))))=[];
[coeff_SI,score_SI,latent_SI,~,explained_SI,mu_SI] = pca(SI_matrix_i);
for l=1:size(Nan_col,2)
    coeff_SI=[coeff_SI([1:Nan_col(1,l)-1],:);zeros(1,size(coeff_SI,2));coeff_SI([Nan_col(1,l):end],:)];
end

for k=1:n_pc
    hf=scatter3(ptCloud_template(:,1),ptCloud_template(:,2),ptCloud_template(:,3),'filled','cdata',coeff_SI(:,k));
    hf=gca; 
    grid off;
    axis off;
    colorbar;
    view([0 90]);
    title([PartName,'/pc',num2str(k)])
    saveas(gca,[Dirpath_output,PartName,'_pc',num2str(k),'.jpeg']);
end
end



