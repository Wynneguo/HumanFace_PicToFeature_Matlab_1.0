function GetAngleFeatures(SampleNames,Path_template_lm,Dirpath_of_Samplelm,Dirpath_output)
    
AngleFeatures={};
template_lm=importdata(Path_template_lm);
lm_mark=template_lm(:,1);%the mark of points

%Rf1 眼睛大小/眼间距
Rf1={};
lm_withindex=[];
for i=1:size(SampleNames,1)
    lm=importdata([Dirpath_of_Samplelm,SampleNames{i,:},'.txt']);
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
    lm=importdata([Dirpath_of_Samplelm,SampleNames{i,:},'.txt']);
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
    lm=importdata([Dirpath_of_Samplelm,SampleNames{i,:},'.txt']);
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
    lm=importdata([Dirpath_of_Samplelm,SampleNames{i,:},'.txt']);
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
    lm=importdata([Dirpath_of_Samplelm,SampleNames{i,:},'.txt']);
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
    lm=importdata([Dirpath_of_Samplelm,SampleNames{i,:},'.txt']);
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
    lm=importdata([Dirpath_of_Samplelm,SampleNames{i,:},'.txt']);
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
    lm=importdata([Dirpath_of_Samplelm,SampleNames{i,:},'.txt']);
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
    lm=importdata([Dirpath_of_Samplelm,SampleNames{i,:},'.txt']);
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
writetable(AngleFeatures,[Dirpath_output,'AngleFeatures.txt'],'Delimiter',' ','WriteVariableNames',true);
%a=readtable([Path_of_features,'features_distance.txt'],'Format',['%s',repmat('%f',1,351)]);
end