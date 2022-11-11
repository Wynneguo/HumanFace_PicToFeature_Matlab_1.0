function [idxSet]=prediction_knn(trainData,trainClass,testData,K)
%% trainData n*m ; trainClass n*1 ; testData s*m
%1）算距离：给定测试对象，计算它与训练集中的每个对象的距离
%2）找邻居：圈定距离最近的k个训练对象，作为测试对象的近邻
%3）做分类：根据这k个近邻归属的主要类别，来对测试对象分类

[N,M]=size(trainData);
S=size(testData,1);
%计算训练数据集与测试数据之间的欧氏距离dist
dist=zeros(N,S);
for j=1:S  
    for i=1:N
        dist(i,j)=norm(trainData(i,:)-testData(j,:));
    end
end
%将dist从小到大进行排序
idxSet=zeros(size(testData,1),1);
for i=1:size(testData,1)
    [Y,I]=sort(dist(:,i),1);   
    K=min(K,length(Y));
%将训练数据对应的类别与训练数据排序结果对应
    labels=trainClass(I);
%{
%确定前K个点所在类别的出现频率
classNum=length(unique(trainClass));%取集合中的单值元素的个数
labels=zeros(1,classNum);
for i=1:K
    j=trainClass(i);
    labels(j)=labels(j)+1;
end
%返回前K个点中出现频率最高的类别作为测试数据的预测分类
[~,idx]=max(labels);
%}
%确定前K个点所在类别的出现频率
    idx=mode(labels(1:K));%mode函数求众数
    idxSet(i)=idx;
end
%fprintf('该测试数据属于类 %d  ',idx);

%% draw
%{
for i=1:size(testData,1)
    if idxSet(i,1)==1
        colorlist=[1 0 1];
    elseif idxSet(i,1)==2
        colorlist=[0 1 1];
    elseif idxSet(i,1)==3
        colorlist=[1 0 0];  
    end
    scatter(testData(i,1),testData(i,2),'MarkerEdgeColor',colorlist);
    hold on;
end
%}
result1=testData(find(idxSet==1),:);
result2=testData(find(idxSet==2),:);
result3=testData(find(idxSet==3),:);
scatter(result1(:,1),result1(:,2),'MarkerEdgeColor',[1 0 1]);
hold on;
scatter(result2(:,1),result2(:,2),'MarkerEdgeColor',[0 1 1]);
hold on;
scatter(result3(:,1),result3(:,2),'MarkerEdgeColor',[1 0 0]);
hold on;
end
