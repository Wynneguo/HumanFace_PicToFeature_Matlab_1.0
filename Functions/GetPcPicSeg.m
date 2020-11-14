function GetPcPicSeg(Dirpath_seg_indice,Path_wholefaceSImatrix,Path_template,n_pc,Dirpath_output)
ptCloud_template=readOBJ(Path_template);

SI_matrix=importdata(Path_wholefaceSImatrix);
[SI_matrix_s,MU_s,SIGMA_s]=zscore(SI_matrix);

File=dir(fullfile(Dirpath_seg_indice,'*.txt'));
IndiceNames={File.name}';
IndiceNames=split(IndiceNames,'.');
IndiceNames=IndiceNames(:,1);
for i=1:size(IndiceNames,1)
    color_data=zeros(size(ptCloud_template,1),n_pc+1);
    color_data(:,1)=1:size(ptCloud_template,1);
    InNames=IndiceNames{i,:};
    indices=importdata([Dirpath_seg_indice,InNames,'.txt']);
    SI_matrix_i=SI_matrix_s(:,indices);
    Nan_col=find(sum(isnan(SI_matrix_i)));
    SI_matrix_i(:,find(sum(isnan(SI_matrix_i))))=[];
    [coeff_SI,score_SI,latent_SI,~,explained_SI,mu_SI] = pca(SI_matrix_i);
    for l=1:size(Nan_col,2)
        coeff_SI=[coeff_SI([1:Nan_col(1,l)-1],:);zeros(1,size(coeff_SI,2));coeff_SI([Nan_col(1,l):end],:)];
    end
    for j=1:size(indices,2)
        color_data(find(color_data(:,1)==indices(1,j)),[2:n_pc+1])=coeff_SI(j,[1:n_pc]);
    end
    for k=1:n_pc
        hf=scatter3(ptCloud_template(:,1),ptCloud_template(:,2),ptCloud_template(:,3),'filled','cdata',color_data(:,k+1));
        hf=gca; 
        grid off;
        axis off;
        colorbar;
        view([0 90]);
        title([InNames,'/pc',num2str(k)])
        saveas(gca,[Dirpath_output,InNames,'_pc',num2str(k),'.jpeg']);
    end
end
        
