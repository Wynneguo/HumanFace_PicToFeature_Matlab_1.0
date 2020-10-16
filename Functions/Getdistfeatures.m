function Getdistfeatures(SampleNames,Path_template_lm,Path_of_Samplelm,Path_output)
%Path_template_lm (27*4)
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
writetable(sample_dist,[Path_output,'features_distance.txt'],'Delimiter',' ','WriteVariableNames',true);
end

